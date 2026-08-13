"""X3: the three permutation nulls, against injected truth.

Each null exists to replace a derivation, so each test injects the thing the derivation assumes
away and checks the null sees it -- and, just as important, that it sees *nothing* when there is
nothing there. A null that fires on clean data is worse than no null.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from permutation_nulls import (  # noqa: E402
    collision,
    column_shuffle,
    curveball,
    d1_pairs,
    graph_nulls,
    hypergeom_sf,
    independence_null,
    linkage_score,
    minor_matrix,
    solve_epsilon,
)

BASES = "ACGT"


def random_barcodes(n, length, rng):
    return sorted({"".join(rng.choice(BASES) for _ in range(length)) for _ in range(n)})


# ----------------------------------------------------------------- A. independence


def test_independent_positions_show_no_excess():
    rng = random.Random(1)
    bcs = random_barcodes(20_000, 8, rng)
    out = independence_null(bcs, 8)
    assert out["predicted_excess_full_length"] == pytest.approx(1.0, abs=0.1)


def test_a_duplicated_position_is_detected():
    """Position 1 copies position 0: the barcode is worth 7 nt, not 8, and every marginal is
    still uniform -- so the marginal product cannot see it and the null must."""
    rng = random.Random(2)
    bcs = sorted({(lambda s: s[0] + s[0] + s[2:])("".join(rng.choice(BASES) for _ in range(8)))
                  for _ in range(20_000)})
    out = independence_null(bcs, 8)
    assert out["predicted_excess_full_length"] > 2.0
    assert out["p_coll_dependent"] > out["p_coll_independent"]


def test_collision_of_a_uniform_position_is_one_quarter():
    import collections

    assert collision(collections.Counter("ACGT" * 250)) == pytest.approx(0.25)


def test_column_shuffle_preserves_every_marginal():
    import collections

    rng = random.Random(3)
    bcs = random_barcodes(5_000, 6, rng)
    before = [collections.Counter(b[j] for b in bcs) for j in range(6)]
    after = column_shuffle(bcs, rng)
    # The shuffle can collide, so compare frequencies rather than counts.
    for j in range(6):
        a = collections.Counter(b[j] for b in after)
        for base in BASES:
            assert a[base] / len(after) == pytest.approx(
                before[j][base] / len(bcs), abs=0.02
            )


# ------------------------------------------------------------- B. the distance-1 graph


def test_d1_pairs_counts_each_pair_once():
    assert d1_pairs(["AAAA", "AAAC", "AAGC", "TTTT"]) == [(0, 1), (1, 2)]
    assert d1_pairs(["AAAA", "AACC"]) == []


def test_solve_epsilon_inverts_its_own_model():
    import math

    counts = [10] * 500
    eps, L = 2e-3, 12
    excess = 3 * L * sum(
        (1 - math.exp(-c * eps / 3)) + (1 - math.exp(-c * eps / 3)) ** 2 for c in counts
    )
    assert solve_epsilon(counts, excess, L) == pytest.approx(eps, rel=0.02)


def test_count_shuffle_finds_injected_parents_and_children():
    """Big parents with small children at distance 1. Nothing about the barcode composition or
    the count distribution changes under the shuffle -- only which count sits on which node."""
    rng = random.Random(4)
    parents = random_barcodes(3_000, 10, rng)[:2_000]
    counts = {p: 60 for p in parents}
    for p in parents[:600]:
        j = rng.randrange(10)
        child = p[:j] + rng.choice([b for b in BASES if b != p[j]]) + p[j + 1:]
        if child not in counts:
            counts[child] = 2
    bcs = sorted(counts)
    out = graph_nulls(bcs, counts, 10, rng, 5, (10,))
    row = out["count_ratio"][0]
    assert row["observed"] > 400
    assert row["excess"] > 0.5 * row["observed"]
    assert row["z"] > 5


def test_count_shuffle_finds_nothing_when_counts_are_unrelated():
    rng = random.Random(5)
    bcs = random_barcodes(20_000, 8, rng)
    counts = {b: rng.randrange(1, 80) for b in bcs}
    out = graph_nulls(bcs, counts, 8, rng, 10, (10,))
    row = out["count_ratio"][0]
    assert abs(row["z"]) < 4


# --------------------------------------------------------------- C. within-MIG linkage


def test_curveball_preserves_both_margins():
    rng = random.Random(6)
    mat = [[1 if rng.random() < 0.3 else 0 for _ in range(7)] for _ in range(20)]
    rows = [sum(r) for r in mat]
    cols = [sum(r[c] for r in mat) for c in range(7)]
    out = curveball(mat, rng, 2_000)
    assert [sum(r) for r in out] == rows
    assert [sum(r[c] for r in out) for c in range(7)] == cols


def test_hypergeom_sf_is_a_survival_function():
    assert hypergeom_sf(0, 10, 4, 4) == pytest.approx(0.0, abs=1e-9)
    assert hypergeom_sf(4, 10, 4, 4) > hypergeom_sf(2, 10, 4, 4) > 0


def test_linked_subclone_scores_far_above_scattered_error():
    """Two molecules in one MIG differ at several positions *on the same reads*. Scattered
    sequencing error puts the same number of minor bases on unrelated reads."""
    rng = random.Random(7)
    ref = "".join(rng.choice(BASES) for _ in range(60))
    alt = list(ref)
    for j in (5, 17, 33, 48):
        alt[j] = next(b for b in BASES if b != ref[j])
    linked = [ref] * 12 + ["".join(alt)] * 8
    scattered = []
    for i in range(20):
        s = list(ref)
        for j in (5, 17, 33, 48):
            if rng.random() < 0.4:
                s[j] = next(b for b in BASES if b != ref[j])
        scattered.append("".join(s))
    s_linked = linkage_score(minor_matrix(linked, 2, 8))
    s_scattered = linkage_score(minor_matrix(scattered, 2, 8))
    assert s_linked > 4.0
    assert s_linked > s_scattered + 2.0


def test_a_bad_read_does_not_look_like_a_subclone_under_the_curveball():
    """Two reads carry a minor base at every callable position -- read-level, not molecule-level.
    The hypergeometric calls that linked; the both-margins null does not, because it keeps those
    reads' error load and only moves *where* it lands."""
    rng = random.Random(8)
    n, k = 24, 6
    mat = [[0] * k for _ in range(n)]
    for i in (0, 1):
        mat[i] = [1] * k
    obs = linkage_score(mat)
    null = [linkage_score(curveball(mat, rng, 400)) for _ in range(30)]
    assert obs > 0
    assert sum(1 for s in null if s >= obs) > 15  # the observed sits mid-null, not in its tail
