"""X2's estimator, checked against a floor we injected ourselves.

The RT/PCR floor is the constant that caps every quality migec emits above ~Q40, so the instrument
that measures it has to be calibrated before it is pointed at real data. Here the simulator injects
a known ``rt_error`` -- present in every read of a molecule, therefore un-removable by any
consensus -- and the estimator has to find it.

This caught the first version of the estimator, which fitted ``p_floor + a/c`` by least squares and
returned a *negative* probability: the 2-read bin, where a majority vote is a coin flip, dominated
the intercept.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from tests.conftest import requires_core
from tests.synthetic._sim import SimConfig, simulate

pytestmark = requires_core

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "quality_floor.py"
_spec = importlib.util.spec_from_file_location("quality_floor", SCRIPT)
qf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(qf)

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def _run(tmp_path, rt_error, seed):
    cfg = SimConfig(
        n_molecules=6000,
        n_clones=1,  # a clonal control: every molecule is the same sequence
        seq_len=140,
        umi_len=12,
        coverage=12.0,
        seq_error=2e-3,
        rt_error=rt_error,
        pcr_error=0.0,
        umi_error=0.0,
        adapter=ADAPTER,
        seed=seed,
    )
    sim = simulate(cfg, tmp_path / "sim")
    out = tmp_path / "x2"
    qf.main(["--reads", str(sim["reads"]), "--out", str(out), "--window", "120"])
    fit = dict(
        line.split("\t") for line in (out / "fit.tsv").read_text().strip().split("\n")
    )
    return {k: float(v) for k, v in fit.items()}


def test_the_injected_floor_is_recovered(tmp_path):
    fit = _run(tmp_path, 1e-4, seed=5)
    assert fit["p_floor_lo"] <= 1e-4 <= fit["p_floor_hi"], fit
    assert fit["p_floor"] == pytest.approx(1e-4, rel=0.5), fit


def test_a_ten_times_lower_floor_is_told_apart_from_the_first(tmp_path):
    # The whole point is settling 1e-4 vs 1e-5 vs 1e-6, so the estimator has to resolve a decade.
    fit = _run(tmp_path, 1e-5, seed=6)
    assert fit["p_floor_lo"] <= 1e-5 <= fit["p_floor_hi"], fit
    assert fit["p_floor_hi"] < 1e-4, f"cannot tell 1e-5 from 1e-4: {fit}"


def test_no_floor_gives_a_bound_rather_than_a_number(tmp_path):
    # With rt_error = 0 the consensus is exact, and the honest output is an upper bound set by how
    # many bases were examined -- not a point estimate of zero dressed up as a measurement.
    fit = _run(tmp_path, 0.0, seed=7)
    assert fit["p_floor"] == 0.0, fit
    assert 0 < fit["p_floor_hi"] < 1e-5, fit
    assert fit["q_cap"] > 45, fit


def test_a_tie_is_not_a_call(tmp_path):
    # Two reads disagreeing have no majority. Resolving that by insertion order turns a coin flip
    # into a confident wrong base, which is exactly what broke the first estimator.
    assert qf.consensus(["ACGT", "ATGT"]) == "ANGT"
    assert qf.consensus(["ACGT", "ATGT", "ATGT"]) == "ATGT"
    assert qf.consensus(["ACGTA", "ACGT"]) == "ACGT"  # shortest common length, no indel guessing


def test_the_interval_covers_zero_counts():
    # A bin with no mismatches still has an upper bound, and it has to shrink as bases accumulate.
    lo_small, hi_small = qf.poisson_ci(0, 10_000)
    lo_big, hi_big = qf.poisson_ci(0, 1_000_000)
    assert lo_small == lo_big == 0.0
    assert hi_big < hi_small
    lo, hi = qf.poisson_ci(50, 1_000_000)
    assert lo < 50 / 1_000_000 < hi
