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


def _peak_rss_in_a_fresh_process(path, out, bucket_bits):
    """peak_rss_bytes() is a process high-water mark, so two runs in one interpreter cannot be
    compared -- the second inherits the first's peak. Each configuration gets its own process."""
    import json
    import subprocess
    import sys

    code = (
        "import json;from migec import _core;"
        f"s=_core.assemble({str(path)!r},{str(out)!r},'S1',1e-4,9.61,False,1,1,{bucket_bits});"
        "print(json.dumps({'rss':s['peak_rss_bytes'],'groups':s['groups'],"
        "'molecules':s['molecules']}))"
    )
    return json.loads(subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                                     check=True).stdout)


def test_peak_memory_is_the_bucket_budget_not_the_library(corpus, tmp_path):
    """The claim that justifies the whole partition: nothing scales with the input.

    One bucket holds the library; sixteen hold a sixteenth. If raising the bucket count does not
    lower peak RSS, the partition is decorative and a real run will be killed by the OOM killer.
    """
    d, path = corpus
    one = _peak_rss_in_a_fresh_process(path, tmp_path / "mem1", 0)
    many = _peak_rss_in_a_fresh_process(path, tmp_path / "mem16", 4)
    print(f"\n   1 bucket:  {one['rss'] / 2**20:6.0f} MB"
          f"\n  16 buckets: {many['rss'] / 2**20:6.0f} MB "
          f"({one['rss'] / max(many['rss'], 1):.2f}x less)")
    assert many["rss"] < one["rss"], "more buckets did not lower peak RSS"
    # The resident part is the bucket, and it is roughly a sixteenth. The floor is the interpreter
    # plus the extension, which is why this is not asserted at 16x.
    assert many["rss"] < 0.75 * one["rss"]
    # ...and the answer is identical either way.
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
