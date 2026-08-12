"""How well checkout groups reads by molecule, scored against the simulator's truth.

This is the migec side of the Calib comparison in ``scripts/compare_calib.py``: the same adjusted
Rand index, on simulated data, so the number is checked on every run rather than only when someone
has Calib installed. Calib clusters on barcode *and* sequence; migec today groups on the barcode
alone, so the residual error here is exactly the UMI collision rate -- which is what makes the
comparison worth running.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from tests.conftest import requires_core
from tests.synthetic._sim import SimConfig, simulate

pytestmark = requires_core

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "compare_calib.py"
_spec = importlib.util.spec_from_file_location("compare_calib", SCRIPT)
cc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cc)

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def _run(tmp_path, cfg):
    from migec.checkout import run

    sim = simulate(cfg, tmp_path / "sim")
    (tmp_path / "barcodes.txt").write_text(f"S1\t{sim['pattern']}\n")
    run(sim["reads"], tmp_path / "barcodes.txt", tmp_path / "out")
    truth = cc.read_truth(Path(sim["truth_reads"]))
    pred = cc.read_migec([tmp_path / "out" / "S1.fq.gz"])
    return sim, cc.adjusted_rand(truth, pred)


def test_grouping_is_near_perfect_on_a_clean_library(tmp_path):
    cfg = SimConfig(
        n_molecules=2000, umi_len=12, coverage=8.0, umi_error=0.0, adapter=ADAPTER, seed=1
    )
    sim, m = _run(tmp_path, cfg)
    # No UMI errors and a 12 nt barcode over 2000 molecules: 4^12 = 16.8 M, so collisions are rare
    # and barcode-only grouping is essentially exact.
    assert m["ari"] > 0.999, m
    assert m["reads_in_split_molecules"] == 0.0, "a clean UMI must never split a molecule"


def test_umi_errors_split_molecules_until_correction_lands(tmp_path):
    # A substitution in the UMI moves a read into a barcode of its own, which is a *split*. This is
    # the failure `migec refine` exists to fix (M3), and Calib avoids it by clustering barcodes at
    # an edit distance. Asserting the direction keeps that honest: the error must show up as
    # splitting, never as merging.
    cfg = SimConfig(
        n_molecules=2000, umi_len=12, coverage=8.0, umi_error=5e-3, adapter=ADAPTER, seed=2
    )
    sim, m = _run(tmp_path, cfg)
    assert m["reads_in_split_molecules"] > 0.0
    assert m["predicted_clusters"] > m["true_molecules"]
    # Merging is what destroys data, and a UMI error cannot cause it.
    assert m["reads_in_merged_clusters"] < 0.02, m


def test_a_short_umi_merges_molecules_and_the_metric_says_so(tmp_path):
    # 6 nt is 4096 barcodes for 2000 molecules -- the birthday bound guarantees collisions, and a
    # collision is a merge. This is the regime where Calib's sequence-aware clustering wins and
    # barcode-only grouping cannot, which is worth having as a measured fact rather than a claim.
    cfg = SimConfig(
        n_molecules=2000, umi_len=6, coverage=8.0, umi_error=0.0, adapter=ADAPTER, seed=3
    )
    sim, m = _run(tmp_path, cfg)
    assert sim["n_umi_collisions"] > 0
    assert m["reads_in_merged_clusters"] > 0.1, m
    assert m["predicted_clusters"] < m["true_molecules"]


def test_the_metric_agrees_with_the_simulator_on_collisions(tmp_path):
    cfg = SimConfig(
        n_molecules=1500, umi_len=8, coverage=6.0, umi_error=0.0, adapter=ADAPTER, seed=4
    )
    sim, m = _run(tmp_path, cfg)
    # Every distinct observed UMI is one predicted cluster, so the two counts must agree exactly.
    assert m["predicted_clusters"] == pytest.approx(sim["n_distinct_umis"], abs=1)
