"""Speed and memory regressions for assemble. Off unless RUN_BENCHMARK=1.

assemble has one memory claim and it is the whole reason it partitions: **one bucket is resident
at a time**, so peak RSS is set by the bucket budget rather than by the library. That is asserted
here against a corpus deliberately built to be memory-hostile -- shallow, so the distinct barcode
count is as large as the read count allows.

Thresholds are loose on purpose. They exist to catch a 10x regression, not a 10% one, because CI
runners vary by more than that.

    RUN_BENCHMARK=1 python -m pytest tests/benchmark -q -s
"""

from __future__ import annotations

import gzip
import os
import random

import pytest

from tests.conftest import requires_core

pytestmark = [
    requires_core,
    pytest.mark.skipif(
        os.environ.get("RUN_BENCHMARK") != "1", reason="set RUN_BENCHMARK=1 to run benchmarks"
    ),
]

N_READS = int(os.environ.get("BENCHMARK_READS", 500_000))
READS_PER_UMI = 4
PAYLOAD = 90


@pytest.fixture(scope="module")
def corpus(tmp_path_factory):
    """A checkout-shaped FASTQ: payload only, barcode in the RX tag."""
    d = tmp_path_factory.mktemp("bench_assemble")
    path = d / "S1.fq.gz"
    rng = random.Random(0)
    with gzip.open(path, "wt", compresslevel=1) as fh:
        i = 0
        while i < N_READS:
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            payload = "".join(rng.choice("ACGT") for _ in range(PAYLOAD))
            for _ in range(READS_PER_UMI):
                fh.write(
                    f"@r{i} RX:Z:{umi}\tQX:Z:{'I' * 12}\tBC:Z:S1\n"
                    f"{payload}\n+\n{'I' * PAYLOAD}\n"
                )
                i += 1
    return d, path


def test_throughput(corpus):
    from migec.assemble import run

    d, path = corpus
    s = run(path, d / "out", gzip_level=1)
    rate = s["reads"] / s["wall_seconds"]
    print(
        f"\n  {rate:,.0f} reads/s = {s['reads'] / s['partition_seconds']:,.0f} partitioning + "
        f"{s['groups'] / (s['wall_seconds'] - s['partition_seconds']):,.0f} groups/s consensus"
    )
    # The consensus is one pass over each column with a table lookup per read base, plus three
    # exps per output base. A tenth of the measured rate means a transcendental moved into the
    # per-read loop or the group is being copied per column.
    assert rate > 20_000, f"assemble throughput collapsed to {rate:,.0f} reads/s"


def _peak_rss_in_a_fresh_process(path, out, bucket_bits, threads=1):
    """peak_rss_bytes() is a process high-water mark, so two runs in one interpreter cannot be
    compared -- the second inherits the first's peak. Each configuration gets its own process."""
    import json
    import subprocess
    import sys

    code = (
        "import json;from migec import _core;"
        f"s=_core.assemble({str(path)!r},{str(out)!r},'S1',1e-4,9.61,False,False,1,1,"
        f"{bucket_bits},{threads});"
        "print(json.dumps({'rss':s['peak_rss_bytes'],'groups':s['groups'],'reads':s['reads'],"
        "'molecules':s['molecules'],'buckets':s['buckets']}))"
    )
    return json.loads(subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                                     check=True).stdout)


def _four_times(path, tmp_path):
    """The same corpus, four times over, so the library grows and the partition does not."""
    import gzip
    import shutil

    out = tmp_path / "four.fq.gz"
    if not out.exists():
        with gzip.open(out, "wb", compresslevel=1) as dst:
            for _ in range(4):
                with gzip.open(path, "rb") as src:
                    shutil.copyfileobj(src, dst)
    return out


def test_peak_memory_is_the_partition_not_the_library(corpus, tmp_path):
    """The claim that justifies the whole partition: RSS does not follow the input size.

    Note: what is resident is one bucket per WORKER, not one bucket. Pass 2 used to be serial, so
    the old form of this test compared bucket counts -- but `bucket_bits=0` has always meant
    "choose", so with the partition now floored at 16 that comparison was 16 buckets against 16
    and asserted nothing. Four times the reads at a fixed thread count is the claim, and it is the
    one a run gets killed for breaking.
    """
    d, path = corpus
    thin = _peak_rss_in_a_fresh_process(path, tmp_path / "thin", 4, threads=1)
    fat = _peak_rss_in_a_fresh_process(_four_times(path, tmp_path), tmp_path / "fat", 4, threads=1)
    print(f"\n  {thin['reads']:>9,} reads: {thin['rss'] / 2**20:6.0f} MB"
          f"\n  {fat['reads']:>9,} reads: {fat['rss'] / 2**20:6.0f} MB "
          f"({fat['rss'] / max(thin['rss'], 1):.2f}x for 4x the reads)")
    assert fat["reads"] == 4 * thin["reads"]
    # 4x the library, nothing like 4x the memory. Loose on purpose: it is here to catch the
    # library being held, not to police bytes.
    assert fat["rss"] < 2.0 * thin["rss"]


def test_threads_cost_memory_and_the_run_says_how_much(corpus, tmp_path):
    """One bucket per worker, so RSS rises with -t. It is bounded and it is reported; a user who
    cannot afford it lowers -t, which is exactly the knob that cannot change the answer."""
    d, path = corpus
    one = _peak_rss_in_a_fresh_process(path, tmp_path / "t1", 4, threads=1)
    many = _peak_rss_in_a_fresh_process(path, tmp_path / "t8", 4, threads=8)
    print(f"\n  16 buckets, 1 thread: {one['rss'] / 2**20:6.0f} MB"
          f"\n  16 buckets, 8 threads: {many['rss'] / 2**20:6.0f} MB")
    assert many["rss"] >= one["rss"]
    # Bounded: eight workers do not cost eight times one, because the bucket is not all of it.
    assert many["rss"] < 8 * one["rss"]
    assert (one["groups"], one["molecules"]) == (many["groups"], many["molecules"])


def test_more_buckets_never_changes_the_output(corpus):
    """The partition is an implementation detail and must not be visible in the result."""
    from migec import _core

    d, path = corpus
    digests = []
    for bits in (0, 2, 5):
        out = d / f"det{bits}"
        _core.assemble(str(path), str(out), "S1", 1e-4, 9.61, False, 1, 1, bits)
        digests.append((out / "S1.mig.tsv").read_bytes())
    assert digests[0] == digests[1] == digests[2]


def test_the_reported_clock_covers_the_whole_run(corpus):
    from migec.assemble import run

    d, path = corpus
    s = run(path, d / "out_clock", gzip_level=1)
    assert 0 < s["partition_seconds"] < s["wall_seconds"]


def test_contig_mode_costs_what_it_costs(corpus):
    """Placement is O(reads^2) inside a group, which is fine at four reads per barcode and is not
    fine at four hundred. The ceiling is asserted so that it is a known one."""
    from migec.assemble import run

    d, path = corpus
    plain = run(path, d / "out_plain", gzip_level=1)
    contig = run(path, d / "out_contig", contig=True, gzip_level=1)
    ratio = contig["wall_seconds"] / plain["wall_seconds"]
    print(f"\n  contig mode is {ratio:.2f}x amplicon at {READS_PER_UMI} reads per barcode")
    assert ratio < 10.0, f"contig placement cost {ratio:.1f}x, not the small multiple expected"


@pytest.fixture(scope="module")
def shallow_corpus(tmp_path_factory):
    """1-3 reads per UMI, the memory-hostile case.

    Distinct barcodes are what everything in assemble scales with -- the sort, the bucket, the
    group loop -- and a shallow library maximises them for a given read count. Bulk repertoire
    profiling and shallow 3' GEX both look like this, so it is the shape to benchmark, not the
    deeply-sequenced amplicon where four reads collapse to one group.
    """
    d = tmp_path_factory.mktemp("bench_shallow")
    path = d / "S1.fq.gz"
    rng = random.Random(1)
    with gzip.open(path, "wt", compresslevel=1) as fh:
        for i in range(N_READS):
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            payload = "".join(rng.choice("ACGT") for _ in range(PAYLOAD))
            fh.write(
                f"@r{i} RX:Z:{umi}\tQX:Z:{'I' * 12}\tBC:Z:S1\n{payload}\n+\n{'I' * PAYLOAD}\n"
            )
    return d, path


def test_shallow_throughput(shallow_corpus):
    from migec.assemble import run

    d, path = shallow_corpus
    s = run(path, d / "out", gzip_level=1)
    rate = s["reads"] / s["wall_seconds"]
    print(f"\n  shallow: {rate:,.0f} reads/s, {s['groups']:,} groups from {s['reads']:,} reads "
          f"({s['reads'] / s['groups']:.2f} reads/UMI)")
    assert s["groups"] > 0.9 * s["reads"], "this fixture is meant to be one read per UMI"
    # Every read is its own group, so the per-group work runs N_READS times instead of N_READS/4.
    # That is the worst case and it still has to clear the same bar.
    assert rate > 20_000, f"shallow throughput collapsed to {rate:,.0f} reads/s"


def test_shallow_memory_is_still_bounded_by_the_bucket(shallow_corpus, tmp_path):
    """The case that would break first if anything held the library: one distinct barcode per read.

    A hash map keyed by barcode would be ~48 B x N_READS here. The partition has to be what bounds
    it, and a finer partition must not cost more than a coarse one at the same thread count.
    """
    d, path = shallow_corpus
    coarse = _peak_rss_in_a_fresh_process(path, tmp_path / "s1", 4, threads=1)
    fine = _peak_rss_in_a_fresh_process(path, tmp_path / "s16", 7, threads=1)
    print(f"\n  shallow  16 buckets, 1 thread: {coarse['rss'] / 2**20:6.0f} MB"
          f"\n  shallow 128 buckets, 1 thread: {fine['rss'] / 2**20:6.0f} MB "
          f"({coarse['rss'] / max(fine['rss'], 1):.2f}x), {fine['groups']:,} groups")
    assert fine["rss"] <= 1.05 * coarse["rss"], "a finer partition cost more memory, not less"
    assert (coarse["groups"], coarse["molecules"]) == (fine["groups"], fine["molecules"])
    many = fine
    # Per distinct barcode, at the finest partition. The sorted (key, count) array in checkout is
    # ~22 B; assemble additionally holds the payload of one bucket, so this is looser -- it is here
    # to catch a hash map reappearing, not to police bytes.
    per_group = many["rss"] / many["groups"]
    print(f"  {per_group:.0f} B per distinct barcode resident")
    assert per_group < 1_000


def test_premade_buckets_skip_the_partition_pass(tmp_path):
    """`checkout --mig` writes the partition, so assemble does not have to.

    The saving is assemble's whole first pass, and it is paid for by a little extra work in
    checkout -- so the number that matters is the two stages together, not either alone. Measured
    on 500,000 reads over four samples at four threads: 1.16 s to 0.98 s for the identical 124,878
    molecules.

    Never: the assertion is that the partition pass is GONE, not that the total is faster. A
    wall-clock comparison across two stages on a shared CI runner is noise; `partition_seconds`
    from a run that did not partition is zero on any machine.
    """
    import gzip

    from migec.assemble import run as assemble
    from migec.checkout import run as checkout

    adapter = "CAGTGGTATCAACGCAGAGT"
    reads = tmp_path / "reads.fq.gz"
    rng = random.Random(0)
    with gzip.open(reads, "wt", compresslevel=1) as fh:
        i = 0
        while i < N_READS:
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            payload = "".join(rng.choice("ACGT") for _ in range(PAYLOAD))
            for _ in range(READS_PER_UMI):
                seq = umi + adapter + payload
                fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
                i += 1
    (tmp_path / "bc.txt").write_text(f"S1\t{'N' * 12}{adapter.lower()}\n")

    fastq = checkout(reads, tmp_path / "bc.txt", tmp_path / "co_fastq", threads=4)
    mig = checkout(reads, tmp_path / "bc.txt", tmp_path / "co_mig", threads=4, mig=True)
    a = assemble(tmp_path / "co_fastq" / "S1.fq.gz", tmp_path / "as_fastq", threads=4)
    b = assemble(tmp_path / "co_mig" / "S1.000.mig", tmp_path / "as_mig", threads=4)
    print(f"\n  FASTQ route: checkout {fastq['wall_seconds']:.2f} s + assemble "
          f"{a['wall_seconds']:.2f} s (partition {a['partition_seconds']:.2f})"
          f"\n  .mig route:  checkout {mig['wall_seconds']:.2f} s + assemble "
          f"{b['wall_seconds']:.2f} s (partition {b['partition_seconds']:.2f})")

    assert b["partition_seconds"] == 0.0, "the partition pass ran on an already-partitioned input"
    assert a["partition_seconds"] > 0.0, "the FASTQ route stopped partitioning -- compare to what?"
    # ...and the point of all of it: the same molecules either way.
    assert b["molecules"] == a["molecules"]
    assert b["groups"] == a["groups"]
