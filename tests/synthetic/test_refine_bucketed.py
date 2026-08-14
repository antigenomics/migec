"""refine's barcode table bounds itself, and bounding it changes nothing but the memory.

Roadmap item 1. `checkout`'s counters already spill; refine's table was the last thing in the
pipeline that scaled with the library, and it is the harder half, because refine's table is not
counts alone -- it carries the per-barcode EVIDENCE (the barcode's own quality at each position
and a draft of its payload) that the posterior needs at 1-3 reads per UMI. A partition that
dropped the evidence would bound the memory and quietly fall back to the count ratio, which
reports nothing in exactly that regime.

So what is asserted here is the whole run against its resident twin: every scalar, and every TSV
byte for byte. A partitioned run that merely *finished* would prove nothing.
"""

from __future__ import annotations

import pytest

from tests.conftest import requires_core
from tests.synthetic._sim import SimConfig, simulate

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"

# Small enough that a few thousand barcodes go over it many times. A table that spilled repeatedly
# is the interesting case: a key is then written to its bucket once per spill, its quality is
# summed across generations, and its payload draft has to stay the FIRST read's however many
# boundaries it crossed.
TINY_BUDGET = 1 << 14

# Every table refine writes. All of them used to be indexed against the entry array; all of them
# are streamed now, so all of them are compared.
TABLES = ["barcodes.tsv", "rank.tsv", "sizes.tsv", "umi_errors.tsv", "bins.tsv"]

# Numbers that describe the clock or the machine rather than the answer.
TIMING = {"wall_seconds", "table_seconds", "correct_seconds", "rewrite_seconds",
          "peak_rss_bytes", "table_resident_bytes", "table_spilled", "threads"}


@pytest.fixture(scope="module")
def library(tmp_path_factory):
    """One simulated library with real barcode errors, refined twice: resident and partitioned."""
    from migec.checkout import run as checkout_run
    from migec.refine import run as refine_run

    d = tmp_path_factory.mktemp("refine_spill")
    cfg = SimConfig(n_molecules=6000, n_clones=30, umi_len=12, coverage=5.0, umi_error=2e-3,
                    adapter=ADAPTER, seed=17)
    sim = simulate(cfg, d / "sim")
    (d / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], d / "bc.txt", d / "co")
    reads = d / "co" / "S1.fq.gz"
    resident = refine_run(reads, d / "resident", table_budget_bytes=0)
    spilled = refine_run(reads, d / "spilled", table_budget_bytes=TINY_BUDGET)
    return d, resident, spilled


def test_the_partition_actually_happened(library):
    _, resident, spilled = library
    assert not resident["table_spilled"]
    assert spilled["table_spilled"]
    # The point of the exercise: what the table HOLDS is the budget, not the library. What it
    # would cost is the same number either way, which is why both are reported.
    assert spilled["table_resident_bytes"] < resident["table_resident_bytes"]
    assert spilled["table_bytes"] == resident["table_bytes"]
    # ...and there was something to correct, or the comparison below is vacuous.
    assert resident["merged"] > 100


def test_every_scalar_matches_the_resident_run(library):
    _, resident, spilled = library
    for key, want in sorted(resident.items()):
        if key in TIMING or not isinstance(want, (int, float, bool)):
            continue
        assert spilled[key] == pytest.approx(want, rel=1e-9), key


def test_the_evidence_survived_the_partition(library):
    """A merge the count ratio alone would have refused is one only the payload could make."""
    _, resident, spilled = library
    assert resident["merged_by_payload"] > 0
    assert spilled["merged_by_payload"] == resident["merged_by_payload"]
    assert spilled["payload_clonality"] == pytest.approx(resident["payload_clonality"], rel=0.25)


@pytest.mark.parametrize("table", TABLES)
def test_every_table_matches_byte_for_byte(library, table):
    d, _, _ = library
    a = (d / "resident" / f"S1.{table}").read_bytes()
    b = (d / "spilled" / f"S1.{table}").read_bytes()
    assert a == b


def test_the_rewritten_reads_match_byte_for_byte(library):
    import gzip

    d, _, _ = library
    a = gzip.decompress((d / "resident" / "S1.fq.gz").read_bytes())
    b = gzip.decompress((d / "spilled" / "S1.fq.gz").read_bytes())
    assert a == b


def test_the_partition_is_removed_after_the_run(library):
    """It is a temporary of the run: leaving it behind would be read back by the next one."""
    d, _, _ = library
    assert not (d / "spilled" / ".refine_spill").exists()
