"""The UMI counters bound themselves, and bounding them changes nothing but the memory.

This is roadmap items 1 and 2, which are one item: the counters range-partition to disk past a
byte budget, and correction -- which walks each barcode's 3L neighbourhood -- follows them into the
partition in two passes with the key rotated. A partition without the rotated pass would bound the
memory and silently stop correcting exactly the children whose substitution landed in the
partitioned prefix, so what is asserted here is the *whole* summary against the resident run, not
just that it ran.
"""

from __future__ import annotations

import pytest

from tests.conftest import requires_core
from tests.synthetic._sim import SimConfig, simulate

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"

# Small enough that a few thousand barcodes go over it many times: the interesting case is a
# counter that spilled repeatedly, because a key is then written to its bucket once per spill and
# reduced only on the way back in.
TINY_BUDGET = 1 << 14


@pytest.fixture(scope="module")
def library(tmp_path_factory):
    """One simulated library with real barcode errors, run twice: resident and partitioned."""
    from migec.checkout import run

    d = tmp_path_factory.mktemp("spill")
    cfg = SimConfig(n_molecules=6000, n_clones=30, umi_len=12, coverage=5.0, umi_error=2e-3,
                    adapter=ADAPTER, seed=17)
    sim = simulate(cfg, d / "sim")
    (d / "bc.txt").write_text(f"S1\t{'N' * cfg.umi_len}{ADAPTER.lower()}\n")
    resident = run(sim["reads"], d / "bc.txt", d / "resident", umi_budget_bytes=0)
    spilled = run(sim["reads"], d / "bc.txt", d / "spilled", umi_budget_bytes=TINY_BUDGET)
    return d, resident, spilled


def test_the_partition_actually_happened(library):
    _, resident, spilled = library
    assert not resident["umi_spilled"]
    assert spilled["umi_spilled"]
    # The point of the exercise: what the counters hold is the budget, not the library. The floor
    # is per counter and there is one sample here, so the resident array is gone and what is left
    # is the append buffer.
    assert spilled["umi_memory_bytes"] < resident["umi_memory_bytes"]


def test_the_partition_is_removed_when_its_readers_are_done(library):
    # The buckets are a temporary of the run: the per-sample statistics stream them *after* the C++
    # returns, so they cannot be deleted at the end of the run -- but a stage that leaves gigabytes
    # next to its output is worse than one that never wrote them.
    d, _, _ = library
    assert list((d / "spilled").rglob(".umi_spill")) == []


def test_every_reported_number_is_the_resident_answer(library):
    _, resident, spilled = library
    r, s = resident["samples"][0], spilled["samples"][0]

    # Counting, composition and the birthday arithmetic stream, so these are exact.
    assert s["umis"] == r["umis"]
    assert s["reads"] == r["reads"]
    assert s["hist_units"] == r["hist_units"]
    assert s["hist_reads"] == r["hist_reads"]
    assert s["effective_length"] == pytest.approx(r["effective_length"])
    assert s["barcode_space"]["molecules"] == pytest.approx(r["barcode_space"]["molecules"])

    # ...and so is correction, which is the claim that costs something to keep true. A one-pass
    # bucketed correction fails here and only here: it merges the children the partition did not
    # cut and misses the rest, which reads as a smaller `umis_merged` and a larger molecule count.
    assert s["umis_merged"] == r["umis_merged"]
    assert s["reads_merged"] == r["reads_merged"]
    assert s["molecules_observed"] == r["molecules_observed"]
    assert s["molecules_corrected"] == pytest.approx(r["molecules_corrected"], rel=1e-6)
    assert s["umi_error_rate"] == pytest.approx(r["umi_error_rate"], rel=0.05)


def test_a_clean_library_reports_the_same_floor(tmp_path):
    """No barcode errors at all: the census finds no excess and there is nothing to solve.

    The resident path then falls back on the 1e-4 floor, so the partitioned one has to report the
    same thing rather than the 0.0 the census returned -- the buckets correct at the floor either
    way, and a summary that says 0.0 is quoting a rate nothing used. It propagates: the error
    budget divides by it.
    """
    from migec.checkout import run

    cfg = SimConfig(n_molecules=5000, n_clones=10, umi_len=12, coverage=4.0, umi_error=0.0,
                    adapter=ADAPTER, seed=3)
    sim = simulate(cfg, tmp_path / "sim")
    (tmp_path / "bc.txt").write_text(f"S1\t{'N' * cfg.umi_len}{ADAPTER.lower()}\n")
    r = run(sim["reads"], tmp_path / "bc.txt", tmp_path / "res", umi_budget_bytes=0)["samples"][0]
    s = run(sim["reads"], tmp_path / "bc.txt", tmp_path / "spl",
            umi_budget_bytes=TINY_BUDGET)["samples"][0]
    assert s["umi_error_rate"] == r["umi_error_rate"]
    assert s["error_budget"]["ratio"] == pytest.approx(r["error_budget"]["ratio"])


def test_the_children_in_the_partitioned_prefix_are_still_corrected(library):
    """The failure mode the rotated pass exists for, stated as a number.

    The partition cuts on the top 8 bits = the first 4 of 12 positions, so a third of the barcode
    errors move a barcode into a different bucket than its parent. If those were lost, `merged`
    would fall by roughly that third -- well outside anything the exact comparison above tolerates.
    """
    _, resident, spilled = library
    r, s = resident["samples"][0], spilled["samples"][0]
    assert r["umis_merged"] > 100, "the library has to contain corrections for this to mean anything"
    assert s["umis_merged"] == r["umis_merged"]
