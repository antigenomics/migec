"""Speed and memory regressions for refine. Off unless RUN_BENCHMARK=1.

refine's memory claim is the opposite of assemble's: it holds the barcode TABLE and streams the
reads, so peak RSS must track *distinct barcodes* and not the read count. The shallow corpus is
the hostile one for exactly that reason -- one distinct barcode per read -- so both are measured.

    RUN_BENCHMARK=1 python -m pytest tests/benchmark -q -s
"""

from __future__ import annotations

import gzip
import json
import os
import random
import subprocess
import sys

import pytest

from tests.conftest import requires_core

pytestmark = [
    requires_core,
    pytest.mark.skipif(
        os.environ.get("RUN_BENCHMARK") != "1", reason="set RUN_BENCHMARK=1 to run benchmarks"
    ),
]

N_READS = int(os.environ.get("BENCHMARK_READS", 500_000))
PAYLOAD = 90
UMI_LEN = 12


def _corpus(path, reads_per_umi, seed, n_barcodes=None):
    """`n_barcodes` fixes the table size independently of the read count, which is what the memory
    claim is actually about."""
    rng = random.Random(seed)
    barcodes = n_barcodes if n_barcodes else N_READS // reads_per_umi
    i = 0
    with gzip.open(path, "wt", compresslevel=1) as fh:
        for _ in range(barcodes):
            umi = "".join(rng.choice("ACGT") for _ in range(UMI_LEN))
            payload = "".join(rng.choice("ACGT") for _ in range(PAYLOAD))
            for _ in range(reads_per_umi):
                fh.write(
                    f"@r{i} RX:Z:{umi}\tQX:Z:{'I' * UMI_LEN}\tBC:Z:S1\n"
                    f"{payload}\n+\n{'I' * PAYLOAD}\n"
                )
                i += 1
    return path


@pytest.fixture(scope="module")
def deep(tmp_path_factory):
    d = tmp_path_factory.mktemp("bench_refine_deep")
    return d, _corpus(d / "S1.fq.gz", 4, 0)


@pytest.fixture(scope="module")
def shallow(tmp_path_factory):
    d = tmp_path_factory.mktemp("bench_refine_shallow")
    return d, _corpus(d / "S1.fq.gz", 1, 1)


def test_throughput(deep):
    from migec.refine import run

    d, path = deep
    s = run(path, d / "out", gzip_level=1)
    rate = s["reads"] / s["wall_seconds"]
    print(f"\n  {rate:,.0f} reads/s over three passes, {s['barcodes']:,} barcodes, "
          f"{s['merged']:,} merged")
    # Three streaming passes plus a 3L neighbourhood walk per barcode. A tenth of the measured
    # rate means the walk stopped being a binary search or a pass was added.
    assert rate > 20_000, f"refine throughput collapsed to {rate:,.0f} reads/s"


def test_the_table_is_what_scales_not_the_reads(deep, shallow):
    """The claim that justifies streaming three times instead of holding the reads."""
    from migec.refine import run

    dd, deep_path = deep
    sd, shallow_path = shallow
    a = run(deep_path, dd / "scale_a", gzip_level=1)
    b = run(shallow_path, sd / "scale_b", gzip_level=1)
    print(f"\n  4 reads/UMI: {a['barcodes']:,} barcodes, table {a['table_bytes'] / 2**20:.0f} MB"
          f"\n  1 read /UMI: {b['barcodes']:,} barcodes, table {b['table_bytes'] / 2**20:.0f} MB")
    # Same read count, 4x the distinct barcodes -> ~4x the table. The reads contribute nothing.
    assert b["barcodes"] > 3 * a["barcodes"]
    assert b["table_bytes"] > 3 * a["table_bytes"]
    # 16 B of (key, count) + a float per barcode position + the payload draft.
    per_barcode = a["table_bytes"] / a["barcodes"]
    print(f"  {per_barcode:.0f} B per distinct barcode")
    assert per_barcode < 16 + 4 * UMI_LEN + 64


def _peak_in_a_fresh_process(path, out, payload_width):
    """peak_rss_bytes is a process high-water mark, so each configuration gets its own process."""
    code = (
        "import json;from migec import _core;"
        f"s=_core.refine({str(path)!r},{str(out)!r},'S1',True,{payload_width > 0},"
        f"{payload_width},0.95,1);"
        "print(json.dumps({'rss':s['peak_rss_bytes'],'table':s['table_bytes'],"
        "'barcodes':s['barcodes'],'reads':s['reads']}))"
    )
    return json.loads(subprocess.run([sys.executable, "-c", code], capture_output=True,
                                     text=True, check=True).stdout)


def test_peak_memory_tracks_barcodes_not_reads(tmp_path):
    """The claim that justifies streaming three times instead of holding the reads.

    Same barcodes, four times the reads. If peak RSS follows the reads, the table is not what is
    being held and the whole design is decoration.
    """
    thin = _corpus(tmp_path / "thin.fq.gz", 2, 3, n_barcodes=100_000)
    fat = _corpus(tmp_path / "fat.fq.gz", 8, 3, n_barcodes=100_000)
    a = _peak_in_a_fresh_process(thin, tmp_path / "a", 32)
    b = _peak_in_a_fresh_process(fat, tmp_path / "b", 32)
    print(f"\n  {a['reads']:,} reads, {a['barcodes']:,} barcodes: RSS {a['rss'] / 2**20:.0f} MB"
          f"\n  {b['reads']:,} reads, {b['barcodes']:,} barcodes: RSS {b['rss'] / 2**20:.0f} MB")
    assert b["reads"] == 4 * a["reads"]
    assert b["barcodes"] == a["barcodes"]
    assert b["table"] == a["table"]
    # Four times the reads, and the resident size barely moves.
    assert b["rss"] < 1.3 * a["rss"]


def test_the_payload_draft_is_a_real_allocation(tmp_path):
    """`table_bytes` has to describe the allocation it claims to, not an arithmetic guess."""
    path = _corpus(tmp_path / "s.fq.gz", 1, 4, n_barcodes=300_000)
    wide = _peak_in_a_fresh_process(path, tmp_path / "w", 64)
    narrow = _peak_in_a_fresh_process(path, tmp_path / "n", 0)
    print(f"\n  payload 64: table {wide['table'] / 2**20:5.0f} MB, RSS {wide['rss'] / 2**20:5.0f} MB"
          f"\n  payload off: table {narrow['table'] / 2**20:5.0f} MB, "
          f"RSS {narrow['rss'] / 2**20:5.0f} MB")
    assert narrow["table"] < wide["table"]
    assert narrow["rss"] < wide["rss"]
    assert wide["table"] - narrow["table"] == 64 * wide["barcodes"]


def test_the_evidence_costs_what_it_costs(deep):
    """Both extra terms need a second pass over the reads and a bigger table. The ceiling is
    asserted so that it is a known one rather than a surprise."""
    from migec.refine import run

    d, path = deep
    full = run(path, d / "cost_full", gzip_level=1)
    counts_only = run(path, d / "cost_counts", use_quality=False, use_payload=False, gzip_level=1)
    ratio = full["wall_seconds"] / counts_only["wall_seconds"]
    print(f"\n  quality + payload cost {ratio:.2f}x wall clock and "
          f"{full['table_bytes'] / max(counts_only['table_bytes'], 1):.1f}x table")
    assert ratio < 3.0, f"the evidence pass cost {ratio:.1f}x, not the small multiple expected"
