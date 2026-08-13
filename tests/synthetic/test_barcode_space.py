"""The birthday arithmetic and the error budget, against injected truth.

These are the two model-derived numbers checkout reports, and both feed decisions: the collision
statistics decide whether a molecule count means anything, the error budget decides whether the
estimated barcode error rate is believable. Each is checked against something the simulator knows.
"""

from __future__ import annotations

import math

import pytest

from tests.conftest import requires_core
from tests.synthetic._sim import SimConfig, simulate

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def _run(tmp_path, cfg, name="out"):
    from migec.checkout import run

    sim = simulate(cfg, tmp_path / f"sim_{name}")
    (tmp_path / f"bc_{name}.txt").write_text(f"S1\t{'N' * cfg.umi_len}{ADAPTER.lower()}\n")
    return run(sim["reads"], tmp_path / f"bc_{name}.txt", tmp_path / name)["samples"][0], sim


def test_a_roomy_barcode_reports_almost_no_collisions(tmp_path):
    cfg = SimConfig(n_molecules=4000, umi_len=12, coverage=6.0, umi_error=0.0,
                    adapter=ADAPTER, seed=1)
    s, sim = _run(tmp_path, cfg)
    b = s["barcode_space"]
    assert b["length"] == 12
    assert b["nominal_space"] == 4**12
    assert b["occupancy"] < 0.001
    assert b["p_multi"] < 0.001
    # With collisions negligible, the corrected molecule count is the observed one.
    assert b["molecules"] == pytest.approx(b["observed"], rel=0.001)
    assert not b["saturated"]


def test_a_short_barcode_reports_the_collisions_the_simulator_made(tmp_path):
    # 4^6 = 4096 barcodes for 3000 molecules. The birthday bound is not subtle here, and the
    # simulator counts the collisions independently.
    cfg = SimConfig(n_molecules=3000, umi_len=6, coverage=6.0, umi_error=0.0,
                    adapter=ADAPTER, seed=2)
    s, sim = _run(tmp_path, cfg)
    b = s["barcode_space"]
    assert b["occupancy"] > 0.4
    assert b["p_multi"] > 0.2
    # `hidden` is the count of molecules no barcode reports. The simulator knows the true number.
    assert b["molecules"] > b["observed"]
    assert b["molecules"] == pytest.approx(cfg.n_molecules, rel=0.15), (b, cfg.n_molecules)
    assert sim["n_umi_collisions"] > 0


def test_lambda_and_p_multi_are_the_poisson_they_claim_to_be(tmp_path):
    cfg = SimConfig(n_molecules=3000, umi_len=7, coverage=6.0, umi_error=0.0,
                    adapter=ADAPTER, seed=3)
    s, _ = _run(tmp_path, cfg)
    b = s["barcode_space"]
    # occupied = S(1 - e^-lambda), and P(k>1 | k>=1) follows from the same lambda.
    assert b["observed"] == pytest.approx(
        b["effective_space"] * (1 - math.exp(-b["lambda"])), rel=1e-6
    )
    occupied = 1 - math.exp(-b["lambda"])
    assert b["p_multi"] == pytest.approx(
        (occupied - b["lambda"] * math.exp(-b["lambda"])) / occupied, rel=1e-6
    )
    assert b["molecules"] == pytest.approx(b["effective_space"] * b["lambda"], rel=1e-6)


def test_a_skewed_composition_costs_barcode_space(tmp_path):
    # 4^L is what a perfect synthesiser would give. A real oligo mix does not, and the loss is
    # what makes collisions more frequent than the length suggests.
    even = SimConfig(n_molecules=3000, umi_len=9, coverage=5.0, umi_error=0.0,
                     adapter=ADAPTER, seed=4)
    skew = SimConfig(n_molecules=3000, umi_len=9, coverage=5.0, umi_error=0.0,
                     adapter=ADAPTER, seed=4, umi_base_freqs=(0.40, 0.10, 0.40, 0.10))
    a, _ = _run(tmp_path, even, "even")
    b, _ = _run(tmp_path, skew, "skew")
    assert a["barcode_space"]["bias_loss"] < 0.02
    assert b["barcode_space"]["bias_loss"] > 0.4, b["barcode_space"]
    assert b["barcode_space"]["effective_space"] < a["barcode_space"]["effective_space"]
    assert b["barcode_space"]["effective_length"] < b["barcode_space"]["length"]
    # ...and the skewed library therefore collides more at the same molecule count.
    assert b["barcode_space"]["p_multi"] > a["barcode_space"]["p_multi"]


def test_a_saturated_space_declines_to_estimate(tmp_path):
    cfg = SimConfig(n_molecules=20000, umi_len=5, coverage=4.0, umi_error=0.0,
                    adapter=ADAPTER, seed=5)
    s, _ = _run(tmp_path, cfg)
    b = s["barcode_space"]
    # Every barcode is occupied, so the inversion would report "no collisions" for the most
    # collided library there can be. It declines instead.
    assert b["saturated"]
    assert b["molecules"] == b["observed"]


def test_the_error_budget_predicts_what_the_estimator_finds(tmp_path):
    # A roomy barcode, so the distance-1 estimator is in its working range, and a UMI error rate
    # far above what the reported Phred alone would give -- so the estimate has to track the
    # injected value rather than the prediction.
    cfg = SimConfig(n_molecules=30000, umi_len=12, coverage=6.0, seq_error=1e-3,
                    umi_error=3e-3, adapter=ADAPTER, seed=6)
    s, _ = _run(tmp_path, cfg)
    e = s["error_budget"]
    assert e["from_phred"] > 0 and e["mean_phred"] > 20
    assert e["from_polymerase"] == pytest.approx(1e-5 * 25)
    assert e["predicted"] == pytest.approx(e["from_phred"] + e["from_polymerase"])
    # The estimator has to land within a factor of two of the 3e-3 that was injected.
    assert 1.5e-3 < e["estimated"] < 6e-3, e
    assert not e["estimate_unreliable"]


def test_the_estimator_admits_when_the_neighbourhood_is_full(tmp_path):
    # The distance-1 estimate is the observed pair count minus the coincidence expectation. Once a
    # barcode's 3L neighbours are mostly real barcodes, that is a small difference of two large
    # numbers and the answer collapses towards zero. Saying so beats reporting it.
    cfg = SimConfig(n_molecules=30000, umi_len=8, coverage=4.0, umi_error=3e-3,
                    adapter=ADAPTER, seed=7)
    s, _ = _run(tmp_path, cfg)
    e = s["error_budget"]
    assert e["neighbour_occupancy"] > 0.05
    assert e["estimate_unreliable"]
    assert e["estimated"] < 3e-3, "the collapse is downward, which is why it must be flagged"


def test_the_phred_prediction_is_the_mean_of_the_probabilities(tmp_path):
    """Not 10^(-mean Q/10). The low-Q tail carries nearly all the error and averaging Q hides it."""
    cfg = SimConfig(n_molecules=3000, umi_len=12, coverage=5.0, mean_qual=30,
                    adapter=ADAPTER, seed=8)
    s, _ = _run(tmp_path, cfg)
    e = s["error_budget"]
    by_hand = sum(q["bases"] * 10 ** (-q["phred"] / 10) for q in s["umi_phred"]) / sum(
        q["bases"] for q in s["umi_phred"]
    )
    assert e["from_phred"] == pytest.approx(by_hand, rel=1e-9)
    # This simulator emits one quality value, so the two forms coincide -- which is itself worth
    # asserting, because it means nothing else is being applied on the way through.
    assert e["from_phred"] == pytest.approx(10 ** (-e["mean_phred"] / 10), rel=1e-9)
    # With any spread they do not, and the difference is the whole reason for averaging the
    # probabilities: 10^(-Q/10) is convex, so the low-Q bases dominate and averaging Q first hides
    # them. Half at Q40 and half at Q10 is an error rate of 5%, not the 0.3% "mean Q25" suggests.
    mixed = 0.5 * 10 ** (-40 / 10) + 0.5 * 10 ** (-10 / 10)
    assert mixed > 10 ** (-25 / 10) * 15
