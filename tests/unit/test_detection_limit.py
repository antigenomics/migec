"""The limit-of-detection arithmetic, which is what decides whether a rare-variant assay can work.

Every number here is checkable by hand, and each test pins a claim the docs make. The one that
matters most is the last: below the pre-amplification floor, more input DNA buys nothing, and an
assay design that ignores that spends money on sequencing which cannot help.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent.parent
spec = importlib.util.spec_from_file_location("detection_limit", ROOT / "scripts" / "detection_limit.py")
dl = importlib.util.module_from_spec(spec)
sys.modules["detection_limit"] = dl
spec.loader.exec_module(dl)


def test_poisson_tail_matches_the_closed_form():
    # P(X >= 1) = 1 - e^-lam, the only one with a one-line closed form.
    for lam in (0.1, 1.0, 5.0):
        assert dl.poisson_at_least(lam, 1) == pytest.approx(1 - math.exp(-lam))
    assert dl.poisson_at_least(0.0, 1) == 0.0
    assert dl.poisson_at_least(0.0, 0) == 1.0


def test_the_lod_inverts_the_tail():
    """LOD is defined by P(Poisson(N f) >= k) = confidence, so feeding it back must round-trip."""
    for molecules, k in [(10_000, 3), (1_000, 1), (500_000, 5)]:
        f = dl.limit_of_detection(molecules, k, 0.95)
        assert dl.poisson_at_least(molecules * f, k) == pytest.approx(0.95, abs=1e-3)


def test_lod_improves_as_one_over_molecules():
    """Ten times the molecules detects ten times lower, which is the whole reason input DNA is
    the lever rather than read depth."""
    a = dl.limit_of_detection(10_000, 3, 0.95)
    b = dl.limit_of_detection(100_000, 3, 0.95)
    assert a / b == pytest.approx(10.0, rel=0.01)


def test_background_counts_one_alternative_allele_not_three():
    """Never: the floor is the chance the base is WRONG; a caller asks about ONE alternative, so
    only a third of those errors imitate the variant being tracked. Using the whole floor
    overstates the background threefold."""
    assert dl.background_molecules(30_000, 1e-4) == pytest.approx(1.0)


def test_pooling_sites_is_what_makes_mrd_possible():
    """Tracking a patient's own variant set multiplies the molecules that can carry evidence.
    Thirty sites is thirty times the evidence and thirty times lower a reachable frequency."""
    r = dl.describe(molecules=30_000, sites=30, floor=1e-9, min_support=3, confidence=0.95)
    assert r["pooled_molecules"] == 900_000
    assert r["lod_pooled"] == pytest.approx(r["lod_one_site"] / 30, rel=0.01)


def test_below_the_floor_more_input_cannot_help():
    """The claim the whole page turns on.

    At the RT floor of 1e-4 a single-strand protocol cannot resolve below ~3.3e-5, because at that
    frequency a true variant molecule is as rare as the chemistry's own false ones. Duplex moves
    the floor by orders of magnitude and the molecule count becomes binding again.
    """
    single = dl.describe(molecules=30_000, sites=30, floor=1e-4, min_support=3, confidence=0.95)
    assert single["vaf_equals_background"] == pytest.approx(3.33e-5, rel=0.01)
    # the molecules promise better than the chemistry can deliver -- that is the trap
    assert single["lod_pooled"] < single["vaf_equals_background"]
    assert single["background_pooled"] > 3, "background should swamp a 3-molecule call"

    duplex = dl.describe(molecules=30_000, sites=30, floor=1e-9, min_support=3, confidence=0.95)
    assert duplex["lod_pooled"] > duplex["vaf_equals_background"], "molecules now bind, not chemistry"
    assert duplex["background_pooled"] < 0.01


def test_genome_equivalents_is_the_haploid_mass():
    # 20 ng / 3.3 pg is ~6,060 haploid copies; x2 strands is the ceiling on molecules per site.
    assert 20 * 1000 / dl.PG_PER_HAPLOID_GENOME == pytest.approx(6060.6, rel=1e-3)


def test_the_named_floors_are_the_documented_ones():
    """These are the same names `migec assemble --rt-error` takes, and drifting from them would
    make the two tools disagree about what `rt` means."""
    assert dl.FLOORS["rt"] == 1e-4
    assert dl.FLOORS["medium"] == 1e-5
    assert dl.FLOORS["high"] == 1e-6
    assert dl.FLOORS["duplex"] < dl.FLOORS["high"]
