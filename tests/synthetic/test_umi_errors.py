"""The barcode error rate read off the CHILDREN, against a known injected rate.

`refine` reports two estimates of the same number and they fail in different directions, so the
test is that each one lands where it is claimed to land -- not that either is always right:

* ``estimated_error`` inverts the excess of distance-1 NEIGHBOURS over what independent draws
  would give. A barcode has only ``3L`` neighbours, so the count saturates and the estimate fails
  downward as the space fills.
* ``error_at_depth`` divides the reads carried by a parent's error children by the ``c * L``
  barcode bases that parent had to miscall. Reads have no ceiling, so this one does not saturate --
  but it counts only children correction actually merged, so it is a lower bound wherever
  correction is incomplete.

Both are bounded by the merges that were made, so on a FULL barcode space both go to zero: there
``correct_umis`` refuses to merge, correctly, because a distance-1 neighbour is more likely a real
molecule than a child. The `saturated` flag is what says the answer is a floor, and this file
asserts that too.

The ratios in the sweep are the ones quoted in `src/refine.cpp` and `docs/umi_errors.rst`.
"""

from __future__ import annotations

import math

import pytest

from tests.conftest import requires_core
from tests.synthetic.test_refine import build

pytestmark = requires_core


@pytest.fixture(scope="module")
def diverse(tmp_path_factory):
    """12 nt barcode at 0.2% occupancy and 25 reads a molecule: the case both estimators fit."""
    d = tmp_path_factory.mktemp("umi_errors_diverse")
    build(d, n_molecules=20_000, n_clones=200, coverage=25.0, coverage_cv=0.5, umi_len=12,
          umi_error=1e-3, seed=11)
    from migec.refine import run

    return run(d / "co" / "S1.fq.gz", d / "ref"), d / "ref"


def test_the_children_recover_the_injected_error_rate(diverse):
    """The headline: on a diverse, deep library it lands on the truth, and so does the Phred."""
    s, _ = diverse
    assert s["error_at_depth"] == pytest.approx(1e-3, rel=0.15)
    # Q30 is -10 log10(1e-3), which is the whole point of reporting a Phred beside it.
    assert s["error_phred"] == pytest.approx(30.0, abs=1.0)
    assert s["error_depth"] == 10
    # ...and it agrees with the independent distance-1 estimate, which is what makes either
    # believable. Never: agreement is the check; neither is the reference for the other.
    assert s["error_at_depth"] == pytest.approx(s["estimated_error"], rel=0.25)


def test_the_all_depth_estimate_is_a_lower_bound(diverse):
    """A child whose parent was never sequenced cannot be merged, so it cannot be counted."""
    s, _ = diverse
    assert 0 < s["error_from_children"] <= s["error_at_depth"] * 1.3


def test_the_table_says_what_the_report_says(diverse):
    """The figure is drawn from the table, so the table has to carry the reported number."""
    s, out = diverse
    rows = (out / "S1.umi_errors.tsv").read_text().strip().split("\n")
    header = rows[0].split("\t")
    assert header == [
        "parent_reads", "parents", "child_barcodes", "child_reads", "children_per_parent",
        "reads_per_parent", "neighbours", "error_from_variants", "error_from_reads",
        "phred_from_reads", "estimate",
    ]
    body = [r.split("\t") for r in rows[1:]]
    assert body, "no depths were tabulated"
    # `neighbours` and `estimate` are constant down the column on purpose: the panel draws both as
    # reference lines, and a figure that needs a value not in its own table will drift from it.
    assert {r[6] for r in body} == {f"{3 * s['umi_length']:.1f}"}
    assert {r[10] for r in body} == {f"{s['error_at_depth']:.6e}"}
    # Every row's own reads estimate is the ratio it claims to be.
    for r in body:
        depth, parents, child_reads = int(r[0]), int(r[1]), int(r[3])
        expected = child_reads / parents / (depth * s["umi_length"])
        assert float(r[8]) == pytest.approx(expected, rel=1e-4)
        if float(r[8]) > 0:
            assert float(r[9]) == pytest.approx(-10 * math.log10(float(r[8])), abs=0.01)


def test_a_deeper_row_holds_more_children_than_a_shallow_one(diverse):
    """The model this inverts, as a monotonicity: c*L bases to miscall means more of them at
    larger c. Asserted on aggregates rather than per row, because one depth is a handful of
    parents and that is exactly why the panel draws points and not a line."""
    _, out = diverse
    rows = [r.split("\t") for r in (out / "S1.umi_errors.tsv").read_text().strip().split("\n")[1:]]
    rows = [(int(r[0]), int(r[1]), int(r[3])) for r in rows if int(r[1]) >= 5]
    assert len(rows) > 20, "not enough depths with parents to compare"
    half = len(rows) // 2
    shallow = sum(cr for _, _, cr in rows[:half]) / sum(p for _, p, _ in rows[:half])
    deep = sum(cr for _, _, cr in rows[half:]) / sum(p for _, p, _ in rows[half:])
    assert deep > shallow


@pytest.mark.parametrize(
    "umi_len, coverage, floor, ceiling",
    [
        (12, 25.0, 0.85, 1.15),   # 0.2% occupancy -- measured 0.99
        (10, 25.0, 0.80, 1.15),   # 2.3%           -- measured 0.95
        (9, 40.0, 0.70, 1.10),    # 9.8%           -- measured 0.88
        (8, 40.0, 0.45, 0.90),    # 33%            -- measured 0.62, and already well low
    ],
)
def test_it_degrades_with_occupancy_rather_than_collapsing(tmp_path, umi_len, coverage, floor,
                                                           ceiling):
    """Occupancy is what breaks a distance-1 census, so the bounds widen with it deliberately.
    They are loose: the point is the trend, not a third decimal place."""
    from migec.refine import run

    build(tmp_path, n_molecules=20_000, n_clones=200, coverage=coverage, coverage_cv=0.5,
          umi_len=umi_len, umi_error=1e-3, seed=11)
    s = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "ref")
    assert floor <= s["error_at_depth"] / 1e-3 <= ceiling


def test_a_full_barcode_space_reports_a_floor_and_says_so(tmp_path):
    """Never: the case that must not be read as a measurement. A 6 nt barcode is 4,096 of them, so
    3,000 molecules fill it; correction then refuses to merge -- rightly, because a neighbour there
    is more likely a real molecule -- and BOTH estimates collapse. What has to survive is the flag
    that says so."""
    from migec.refine import run

    build(tmp_path, n_molecules=3_000, n_clones=100, coverage=60.0, coverage_cv=0.6, umi_len=6,
          umi_error=5e-3, seed=3)
    s = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "ref")
    assert s["saturated"], "a full barcode space must be flagged"
    assert s["error_at_depth"] < 5e-3, "a collapsed estimate must not be mistaken for the truth"
    assert s["estimated_error"] < 5e-3
