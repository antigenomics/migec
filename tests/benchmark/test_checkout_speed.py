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
            sample = list(TAGS)[i % 4]
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
    # Anything below 2x on four cores means work migrated back onto the serial path.
    assert speedup > 2.0, f"4 threads gave only {speedup:.2f}x on matching"


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
