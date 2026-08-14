"""Speed and memory regressions for checkout. Off unless RUN_BENCHMARK=1.

These are the two numbers that decide whether a pipeline can be run at all, so they get asserted
rather than eyeballed. The thresholds are deliberately loose -- they exist to catch a 10x
regression, not to police a 10% one, because CI runners vary by more than that.

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

BARCODES = """\
S1\taaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
S2\taaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN
S3\taaGCCcagtggtatcaacgcagagtNNNNtNNNNtNNNN
S4\taaGGTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
"""
TAGS = {"S1": "ACT", "S2": "AGA", "S3": "GCC", "S4": "GGT"}
ADAPTER = "CAGTGGTATCAACGCAGAGT"

N_READS = int(os.environ.get("BENCHMARK_READS", 500_000))
# Reads per molecule. Low on purpose: the UMI counters scale with *distinct* UMIs, so a library
# sequenced shallowly is the memory-hostile case, not the deep one.
READS_PER_UMI = 4


@pytest.fixture(scope="module")
def corpus(tmp_path_factory):
    d = tmp_path_factory.mktemp("bench")
    path = d / "reads.fq.gz"
    rng = random.Random(0)
    with gzip.open(path, "wt", compresslevel=1) as fh:
        i = 0
        while i < N_READS:
            # On the MOLECULE, not on the read counter: `i % 4` with four reads per molecule is
            # always 0, so every read went to S1 and three of the four patterns never matched.
            # The matching cost was right -- all four are still scored per read -- but the
            # per-sample side of it was measuring one sample.
            sample = list(TAGS)[(i // READS_PER_UMI) % 4]
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            payload = "".join(rng.choice("ACGT") for _ in range(90))
            for _ in range(READS_PER_UMI):
                seq = (
                    "AA" + TAGS[sample] + ADAPTER + umi[:4] + "T" + umi[4:8] + "T" + umi[8:]
                    + payload
                )
                fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
                i += 1
    (d / "barcodes.txt").write_text(BARCODES)
    return d, path


def test_single_thread_throughput(corpus):
    from migec.checkout import run

    d, path = corpus
    s = run(path, d / "barcodes.txt", d / "out1", threads=1)
    rate = s["reads_per_second"]
    print(f"\n  1 thread: {rate:,.0f} reads/s end to end, "
          f"{s['total'] / s['match_seconds']:,.0f} matching")
    # Measured ~193k reads/s end to end on an M-series laptop for 129 nt reads over four patterns.
    # A tenth of that means something structural broke -- a transcendental back in the scoring
    # loop, or compression back on the serial path.
    assert rate > 20_000, f"single-thread throughput collapsed to {rate:,.0f} reads/s"


def test_threads_actually_help(corpus):
    from migec.checkout import run

    d, path = corpus
    # On `match_seconds`, not the end-to-end clock: the per-sample statistics are serial by
    # construction, so including them measures Amdahl's law rather than whether the workers work.
    one = run(path, d / "barcodes.txt", d / "out_a", threads=1)
    many = run(path, d / "barcodes.txt", d / "out_b", threads=4)
    speedup = one["match_seconds"] / many["match_seconds"]
    print(f"\n  4 threads: {speedup:.2f}x on matching, "
          f"{one['wall_seconds'] / many['wall_seconds']:.2f}x end to end")
    # Matching and compression both run on the workers, so the serial part is fread plus fwrite.
    #
    # Never: the floor is what a SHARED four-vCPU runner can meet, not what a quiet laptop does.
    # This asked for 2.0x and a nightly run produced 1.93x -- the runner's four vCPUs are
    # hyperthread siblings shared with other tenants, so the fourth worker is not a fourth core.
    # 1.5x still fails the thing this test is for, which is work migrating back onto the serial
    # path (that shows up as ~1.0x), and it stops the nightly flapping on a number about the
    # machine rather than about the code. Measured on a quiet laptop: 3.3x.
    assert speedup > 1.5, f"4 threads gave only {speedup:.2f}x on matching"


def test_the_reported_clock_covers_the_whole_run(corpus):
    """`reads_per_second` must describe checkout, not the part of it that threads.

    The per-sample statistics -- histogram, composition, count correction -- are a serial pass over
    every distinct UMI, and on a shallow library they are most of the wall clock. A stopwatch that
    stopped at the demultiplexing driver would report several times the throughput a user sees.
    """
    from migec.checkout import run

    d, path = corpus
    s = run(path, d / "barcodes.txt", d / "out_clock", threads=4)
    assert s["match_seconds"] < s["wall_seconds"], "the stats stage is outside the clock again"
    assert s["reads_per_second"] == pytest.approx(s["total"] / s["wall_seconds"], rel=1e-6)


def test_output_is_independent_of_thread_count(corpus):
    from migec.checkout import run

    d, path = corpus
    digests = []
    for t in (1, 2, 8):
        out = d / f"det{t}"
        run(path, d / "barcodes.txt", out, threads=t)
        digests.append((out / "S1.fq.gz").read_bytes())
    assert digests[0] == digests[1] == digests[2]


def test_memory_per_distinct_umi(corpus):
    from migec.checkout import run

    d, path = corpus
    s = run(path, d / "barcodes.txt", d / "out_mem", threads=4)
    distinct = sum(x["umis"] for x in s["samples"])
    per_umi = s["umi_memory_bytes"] / distinct
    print(f"\n  {distinct:,} distinct UMIs, {per_umi:.1f} B each, "
          f"peak RSS {s['peak_rss_bytes'] / 2**20:.0f} MB")
    # The sorted (key, count) array is 16 B per entry plus slack from the growth doubling. A hash
    # map of the same contents runs to ~48 B, which is the difference between a NovaSeq run fitting
    # in memory and not.
    assert per_umi < 48, f"{per_umi:.1f} bytes per distinct UMI -- the flat array regressed"


def test_peak_memory_does_not_track_input_size(corpus):
    """Everything but the UMI counters is bounded by chunk x threads, not by the file."""
    from migec.checkout import run

    d, path = corpus
    s = run(path, d / "barcodes.txt", d / "out_bound", threads=4)
    overhead = s["peak_rss_bytes"] - s["umi_memory_bytes"]
    print(f"\n  non-counter RSS {overhead / 2**20:.0f} MB")
    assert overhead < 1 << 30, "buffers should be bounded by chunk size x threads"


def _corpus_of(directory, reads):
    """A shallow library of `reads` reads, one distinct barcode per 4 of them."""
    path = directory / f"reads_{reads}.fq.gz"
    rng = random.Random(0)
    with gzip.open(path, "wt", compresslevel=1) as fh:
        i = 0
        while i < reads:
            sample = list(TAGS)[(i // READS_PER_UMI) % 4]
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            payload = "".join(rng.choice("ACGT") for _ in range(90))
            for _ in range(READS_PER_UMI):
                seq = ("AA" + TAGS[sample] + ADAPTER + umi[:4] + "T" + umi[4:8] + "T" + umi[8:]
                       + payload)
                fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
                i += 1
    (directory / "barcodes.txt").write_text(BARCODES)
    return path


def test_the_counters_do_not_grow_with_the_library(tmp_path):
    """The allocation that used to grow with the library, stated as something that can fail.

    Everything else in `checkout` is bounded by chunk x threads. The UMI counters were not: ~22 B
    per DISTINCT barcode held in one piece, which is 8.8 GB at NovaSeq scale. They range-partition
    to disk past `umi_budget_bytes` now, and what this asserts is the property that fix has and a
    warning did not: doubling the library does not double the counters.

    Note: the budget is passed explicitly and is small. The default is 1 GB, which no corpus that
    fits in CI can reach -- so a run at the default would pass this test resident, vacuously, and
    the assertion would be about nothing. The budget is the axis being tested, not the corpus.

    Never: assert a SCALING property, not a fixed budget. A budget test passes vacuously at any
    corpus small enough to run in CI -- 50,000 barcodes is 1 MB, under any threshold worth naming,
    so the test would go green while the term it exists to bound grew without limit. Doubling the
    distinct barcodes and watching what the counters do is the question that has the same answer at
    every scale.

    Never: the fix is the range partition, not a smaller struct. 16 B per entry is already near the
    floor for (key, count); shaving it buys a constant factor against a term that grows without
    limit, which is the wrong axis.
    """
    from migec.checkout import run

    budget = 1 << 16  # 64 kB over 4 samples, against ~0.4 MB and ~1.6 MB of barcodes
    small = run(_corpus_of(tmp_path, 100_000), tmp_path / "barcodes.txt", tmp_path / "a",
                threads=4, umi_budget_bytes=budget)
    large = run(_corpus_of(tmp_path, 400_000), tmp_path / "barcodes.txt", tmp_path / "b",
                threads=4, umi_budget_bytes=budget)
    assert small["umi_spilled"] and large["umi_spilled"], "the budget was never reached"

    def distinct(s):
        return sum(x["umis"] for x in s["samples"])

    growth = large["umi_memory_bytes"] / max(small["umi_memory_bytes"], 1)
    print(f"\n  {distinct(small):,} barcodes: {small['umi_memory_bytes'] / 2**20:6.1f} MB"
          f" ({small['umi_memory_bytes'] / distinct(small):.1f} B each)"
          f"\n  {distinct(large):,} barcodes: {large['umi_memory_bytes'] / 2**20:6.1f} MB"
          f" ({large['umi_memory_bytes'] / distinct(large):.1f} B each)"
          f"  ({growth:.2f}x for {distinct(large) / distinct(small):.2f}x the barcodes)")
    # Note: not 1.0x. What is left is a constant, not a term in the library: four append buffers of
    # up to 4096 entries each, plus whatever sat in the sorted array when the input ended, and both
    # are quantised by the vector's growth doubling. Measured 1.78x for 4x the barcodes, i.e. 12.6
    # B per barcode falling to 5.2. Never: the assertion is on the SHAPE -- sub-linear at every
    # scale -- because a fixed budget passes vacuously on any corpus that fits in CI.
    assert growth < 2.0, (
        f"4x the distinct barcodes cost {growth:.2f}x the counter memory -- the partition is not "
        f"bounding them")
