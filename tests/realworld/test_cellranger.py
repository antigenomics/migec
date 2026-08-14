"""The Cell Ranger side of `scripts/compare_cellranger.py`, asserted on its own.

Only the comparator is checked here, not migec: the published Cell Ranger 5.0.0 output for
`sc5p_v2_hs_PBMC_1k` VDJ-T is small, committed nowhere, and fetched on demand, while migec's side
needs the 1.1 GB FASTQ pair. What this tier is for is catching the day 10x reissue the files or
change a column's meaning -- which is exactly the failure the script's read-back exists to refuse,
so the read-back itself is what gets tested.

Note: the CI fixture is a 1% cell subsample and CANNOT show cell calling with a knee, so nothing
here asserts a migec cell count. That belongs on the full library.
"""

from __future__ import annotations

import csv
import sys
import urllib.request
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

BASE = ("https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_1k/"
        "sc5p_v2_hs_PBMC_1k_t")
FILES = ("filtered_contig_annotations.csv", "airr_rearrangement.tsv", "metrics_summary.csv")


@pytest.fixture(scope="module")
def published(tmp_path_factory):
    """10x's own Cell Ranger 5.0.0 output, fetched once. Skips when offline."""
    d = tmp_path_factory.mktemp("cellranger")
    for name in FILES:
        # cf.10xgenomics.com answers 403 to urllib's default User-Agent but serves curl's, so the
        # header is required rather than polite -- without it every test here skips and the suite
        # reports green having checked nothing.
        req = urllib.request.Request(f"{BASE}_{name}", headers={"User-Agent": "curl/8"})
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                (d / f"sc5p_v2_hs_PBMC_1k_t_{name}").write_bytes(r.read())
        except OSError as exc:
            pytest.skip(f"cannot reach cf.10xgenomics.com: {exc}")
    return d


def test_the_cell_set_is_the_one_the_comparison_was_scored_against(published):
    from compare_cellranger import CELLRANGER_CELLS, CELLRANGER_CONTIGS, read_cellranger_cells

    cells, contigs = read_cellranger_cells(
        published / "sc5p_v2_hs_PBMC_1k_t_filtered_contig_annotations.csv")
    assert len(cells) == CELLRANGER_CELLS
    assert contigs == CELLRANGER_CONTIGS
    assert all(len(b) == 16 and set(b) <= set("ACGTN") for b in cells)


def test_is_cell_actually_filters(published):
    """`_all_contig_annotations` spans 699 barcodes; reading it unfiltered inflates the set 46%."""
    path = published / "sc5p_v2_hs_PBMC_1k_t_filtered_contig_annotations.csv"
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert {r["is_cell"].strip().lower() for r in rows} == {"true"}


def test_the_airr_count_columns_still_mean_what_the_script_assumes(published):
    """consensus_count is READS and duplicate_count is UMIs -- the inverse of arda's spelling.

    This is the one that would corrupt a join silently, so it is the one with a test.
    """
    from compare_cellranger import check_airr_counts

    check_airr_counts(published / "sc5p_v2_hs_PBMC_1k_t_airr_rearrangement.tsv",
                      published / "sc5p_v2_hs_PBMC_1k_t_filtered_contig_annotations.csv")


def test_the_headline_metrics_are_the_ones_quoted_in_the_docs(published):
    from compare_cellranger import _pct, read_metrics

    m = read_metrics(published / "sc5p_v2_hs_PBMC_1k_t_metrics_summary.csv")
    assert _pct(m["Valid Barcodes"]) == pytest.approx(0.906)
    assert _pct(m["Fraction Reads in Cells"]) == pytest.approx(0.868)
    assert int(m["Number of Read Pairs"].replace(",", "")) == 6_301_573
