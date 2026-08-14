"""`--pre-amp-error auto`: fit the floor from the data, and refuse when the data cannot carry it.

This file is the calibration of that estimator, so the table in `docs/quality_floor.rst` is these
assertions. Two things are checked, and the second matters as much as the first: that an injected
floor comes back inside its interval, and that a library which cannot support the fit -- diverse,
or too shallow to have deep molecules -- is REFUSED with a named fallback rather than answered
with its own diversity.
"""

from __future__ import annotations

import gzip

import pytest

pytest.importorskip("migec._core", reason="the C++ extension is not built: run `bash setup.sh`")

from migec.assemble import format_report, run
from migec.checkout import run as checkout_run

from ._sim import SimConfig, simulate

from tests.conftest import requires_core

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def fit(tmp_path, **kwargs):
    """Simulate, check out, assemble with `auto`. Returns (summary, estimate)."""
    cfg = SimConfig(adapter=ADAPTER, **kwargs)
    sim = simulate(cfg, tmp_path / "sim")
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co")
    summary = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "asm", rt_floor="auto")
    return summary, summary["pre_amp_estimate"]


def test_an_injected_floor_comes_back_inside_its_interval(tmp_path):
    summary, e = fit(
        tmp_path, n_clones=1, n_molecules=2000, coverage=40.0, coverage_cv=0.4,
        rt_error=1e-4, seed=7,
    )
    assert e["ok"], e["reason"]
    assert e["ci_lo"] <= 1e-4 <= e["ci_hi"]
    # The point estimate sits slightly above the injected RT rate because the simulator also puts
    # early-PCR errors into the molecule, and those are pre-amplification error too -- the floor
    # is what survives a consensus, whatever made it.
    assert 1e-4 <= e["rate"] <= 3e-4
    assert summary["rt_floor"] == e["rate"]
    # Every input the fit rested on is on disk next to the number.
    row = (tmp_path / "asm" / "assemble.pre_amp_error.tsv").read_text().splitlines()
    assert len(row) == 2 and str(e["mismatches"]) in row[1]
    # The probe assembly is scratch and does not survive the run.
    assert not (tmp_path / "asm" / ".pre_amp_probe").exists()
    # And the real output is still a whole file.
    with gzip.open(tmp_path / "asm" / "S1.consensus.fq.gz", "rt") as fh:
        assert sum(1 for _ in fh) % 4 == 0
    assert "pre-amp floor fitted" in format_report(summary)


def test_a_decade_lower_is_resolved_as_a_decade_lower(tmp_path):
    _, e = fit(
        tmp_path, n_clones=1, n_molecules=8000, coverage=30.0, coverage_cv=0.4,
        rt_error=1e-5, pcr_error=1e-6, seed=8,
    )
    assert e["ok"], e["reason"]
    assert e["ci_lo"] <= 1e-5 <= e["ci_hi"]
    assert e["rate"] < 1e-4  # the whole point: it does not answer with the default class


def test_no_injected_error_gives_a_bound_rather_than_zero(tmp_path):
    summary, e = fit(
        tmp_path, n_clones=1, n_molecules=4000, coverage=30.0, coverage_cv=0.4,
        rt_error=0.0, pcr_error=0.0, seed=9,
    )
    assert e["ok"] and e["mismatches"] == 0 and e["rate"] == 0.0
    # Never: a floor of zero is a Q-infinity. With no observed error the answer is the upper end
    # of the interval, which is what the data can actually support.
    assert summary["rt_floor"] == e["ci_hi"] > 0.0
    assert summary["rt_floor"] < 2e-5


def test_a_diverse_library_is_refused_rather_than_told_its_own_diversity(tmp_path):
    summary, e = fit(
        tmp_path, n_clones=400, n_molecules=2000, coverage=40.0, coverage_cv=0.4,
        rt_error=1e-4, seed=10,
    )
    assert not e["ok"] and e["fallback"] == "rt"
    assert "clonal" in e["reason"]
    assert summary["rt_floor"] == 1e-4
    assert "NOT fitted" in format_report(summary)


def test_a_shallow_library_is_refused_because_the_curve_never_flattens(tmp_path):
    summary, e = fit(
        tmp_path, n_clones=1, n_molecules=4000, coverage=1.2, coverage_cv=0.8,
        rt_error=1e-4, seed=11,
    )
    assert not e["ok"] and e["molecules_deep"] == 0
    assert "20 reads or more" in e["reason"]
    assert summary["rt_floor"] == 1e-4
