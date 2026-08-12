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


def test_real_variation_is_excluded_rather_than_counted_as_error():
    """A quasispecies position must not be scored as error.

    The clonal-control assumption is that any disagreement with the modal base is error we
    introduced. On an HIV plasma population that is false: a position carrying a real 20% variant
    would contribute 0.2 to the "error" rate and swamp a floor of 1e-4 by three orders of
    magnitude. The floor is measured only where the molecules agree.
    """
    # 2000 molecules: position 1 carries a genuine 20% variant, position 3 one erroneous molecule.
    cons = []
    for i in range(2000):
        cons.append("A" + ("C" if i < 400 else "A") + "A" + ("G" if i == 0 else "A"))
    spectrum = qf.minor_allele_spectrum(cons, 4)
    assert spectrum[1] == pytest.approx(0.20)
    assert spectrum[3] == pytest.approx(1 / 2000)
    assert spectrum[0] == spectrum[2] == 0.0

    mono = [j for j, m in enumerate(spectrum) if m < 0.01]
    assert mono == [0, 2, 3], "the 20% variant is excluded, the 1-molecule error is kept"


def test_the_threshold_must_sit_above_one_molecule(tmp_path):
    """The exclusion threshold cannot be so tight that a single erroneous molecule looks real.

    With M molecules one error is a minor fraction of 1/M. If ``--max-minor`` is at or below that,
    every position where an error occurred is excluded as "real variation" and the floor is
    measured only where nothing went wrong -- biased down, which is the direction that lets migec
    emit a quality it cannot support.
    """
    from tests.synthetic._sim import SimConfig, simulate

    cfg = SimConfig(
        n_molecules=600, n_clones=1, seq_len=140, umi_len=12, coverage=10.0,
        seq_error=2e-3, rt_error=1e-4, pcr_error=0.0, umi_error=0.0, adapter=ADAPTER, seed=9,
    )
    sim = simulate(cfg, tmp_path / "sim")
    with pytest.raises(SystemExit, match="which is the signal itself"):
        qf.main(["--reads", str(sim["reads"]), "--out", str(tmp_path / "x2"),
                 "--window", "120", "--max-minor", "0.001"])


def test_a_tie_does_not_vote_in_the_spectrum():
    # An unresolved consensus base is not evidence for or against a variant.
    assert qf.minor_allele_spectrum(["A", "N", "A"], 1)[0] == 0.0
    assert qf.minor_allele_spectrum(["N", "N"], 1)[0] == 0.0
