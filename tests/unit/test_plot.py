"""`migec plot`: gnuplot scripts over the tables the stages wrote.

The scripts are the deliverable, so they are what is tested. Rendering is checked too, but only
where gnuplot exists -- it is not a Python package and CI need not have one.
"""

from __future__ import annotations

import shutil

import pytest

from migec.plot import PANELS, format_report, run

COMPOSITION = """sample_id\tposition\tA\tC\tG\tT\tentropy_bits\tinformation_bits\tcollision
S1\t0\t0.25\t0.25\t0.25\t0.25\t2.000000\t0.000000\t0.250000
S1\t1\t0.40\t0.20\t0.20\t0.20\t1.921928\t0.078072\t0.280000
S2\t0\t0.10\t0.30\t0.30\t0.30\t1.895462\t0.104538\t0.280000
"""

COVERAGE = """sample_id\tmig_size\treads\tunits
S1\t1\t100\t100
S1\t2\t400\t200
"""


def test_a_script_is_written_for_every_table_that_exists(tmp_path):
    (tmp_path / "checkout.umi_composition.tsv").write_text(COMPOSITION)
    (tmp_path / "checkout.coverage.tsv").write_text(COVERAGE)

    summary = run(tmp_path, gnuplot="")  # scripts only, as on a machine without gnuplot
    scripts = set(summary["scripts"])
    # One figure per sample for the per-sample panels: two samples in the table, two scripts.
    assert "umi_pwm.S1.gp" in scripts
    assert "umi_pwm.S2.gp" in scripts
    assert "coverage.gp" in scripts
    # ...and nothing invented for the tables that are not there.
    assert "consensus_quality.gp" not in scripts
    assert "cell_rank" in summary["skipped"]
    assert summary["drawn"] == []

    text = (tmp_path / "plots" / "umi_pwm.S1.gp").read_text()
    assert 'strcol(1) eq "S1"' in text          # the filter, so S2's rows stay out of S1's figure
    assert str(tmp_path / "checkout.umi_composition.tsv") in text  # absolute: gnuplot chdirs
    assert "gnuplot was not found" in format_report(summary)


def test_pointing_it_at_a_file_says_so(tmp_path):
    f = tmp_path / "checkout.coverage.tsv"
    f.write_text(COVERAGE)
    with pytest.raises(ValueError, match="is not a directory"):
        run(f)


def test_an_unknown_format_names_the_ones_that_work(tmp_path):
    (tmp_path / "checkout.coverage.tsv").write_text(COVERAGE)
    with pytest.raises(ValueError, match="svg"):
        run(tmp_path, fmt="jpeg")


def test_every_panel_has_a_unique_name_and_a_plot_command():
    names = [p.name for p in PANELS]
    assert len(names) == len(set(names))
    for p in PANELS:
        assert "plot " in p.script
        assert "{src}" in p.script


@pytest.mark.skipif(shutil.which("gnuplot") is None, reason="gnuplot is not installed")
def test_gnuplot_actually_renders_them(tmp_path):
    (tmp_path / "checkout.umi_composition.tsv").write_text(COMPOSITION)
    (tmp_path / "checkout.coverage.tsv").write_text(COVERAGE)

    summary = run(tmp_path, fmt="svg")
    assert summary["failed"] == []
    assert set(summary["drawn"]) == {"umi_pwm.S1.svg", "umi_pwm.S2.svg",
                                     "umi_information.S1.svg", "umi_information.S2.svg",
                                     "coverage.svg"}
    for name in summary["drawn"]:
        assert (tmp_path / "plots" / name).stat().st_size > 1000


SIZES = """size\tlog1p_size\tmolecules\treads
1\t0.693147\t500\t500
2\t1.098612\t120\t240
5\t1.791759\t20\t100
"""

CELL_RANK = """rank\tumis\tcalled\tcumulative_umis\tcumulative_fraction
1\t200\t1\t200\t0.400000
2\t150\t1\t350\t0.700000
3\t2\t0\t352\t0.704000
"""

QUALITY_BY_DEPTH = """sample_id\tmin_reads\tmax_reads\tmolecules\tq_min\tq_p25\tq_median\tq_p75\tq_max\tq_mean
S1\t1\t1\t500\t28\t30\t31\t31\t34\t30.760
S1\t2\t3\t120\t37\t39\t40\t40\t40\t39.609
"""


def test_the_familiar_panels_draw_from_their_own_tables(tmp_path):
    """The four figures a user already knows how to read, and the tables they come off."""
    (tmp_path / "S1.sizes.tsv").write_text(SIZES)
    (tmp_path / "S1.cell_rank.tsv").write_text(CELL_RANK)
    (tmp_path / "assemble.quality_by_depth.tsv").write_text(QUALITY_BY_DEPTH)

    scripts = set(run(tmp_path, gnuplot="")["scripts"])
    assert {"mig_size_spectrum.gp", "mig_size_zipf.gp", "cell_rank.gp",
            "consensus_quality.gp"} <= scripts

    # The knee plot is on UMIs, never reads: column 2 of the cell-rank table is `umis`, and the
    # call (column 3) is what splits the two colours. Drawing reads here would hide the exact
    # artefact the plot exists to show.
    knee = (tmp_path / "plots" / "cell_rank.gp").read_text()
    assert "unique UMIs" in knee
    assert "$3 == 1" in knee and "$3 == 0" in knee

    # A box, not a thinned scatter: `every` would mean a sample of the molecules was drawn.
    quality = (tmp_path / "plots" / "consensus_quality.gp").read_text()
    assert "candlesticks" in quality
    # No `every N`: that is a SAMPLE of the molecules drawn as if it were the distribution.
    plot_line = quality[quality.index("plot "):]
    assert " every " not in plot_line


def test_the_figures_are_transparent_and_one_ink_colour(tmp_path):
    """One SVG has to serve a light page, a dark page and print -- so no background, one grey."""
    (tmp_path / "S1.cell_rank.tsv").write_text(CELL_RANK)
    script = (tmp_path / "plots" / "cell_rank.gp")
    run(tmp_path, gnuplot="")
    text = script.read_text()
    assert "background rgb" not in text      # no background rect => transparent SVG
    assert 'set key inside' in text          # never a legend gutter: it widens every figure
    assert text.count("#808080") >= 4        # border, tics, labels, key


@pytest.mark.skipif(shutil.which("gnuplot") is None, reason="gnuplot is not installed")
def test_gnuplot_renders_the_familiar_panels(tmp_path):
    (tmp_path / "S1.sizes.tsv").write_text(SIZES)
    (tmp_path / "S1.cell_rank.tsv").write_text(CELL_RANK)
    (tmp_path / "assemble.quality_by_depth.tsv").write_text(QUALITY_BY_DEPTH)

    summary = run(tmp_path)
    assert summary["failed"] == []
    assert {"mig_size_spectrum.svg", "mig_size_zipf.svg", "cell_rank.svg",
            "consensus_quality.svg"} <= set(summary["drawn"])
    svg = (tmp_path / "plots" / "cell_rank.svg").read_text()
    assert 'width="760"' in svg
    assert 'fill="none"' in svg              # the canvas rect, unfilled
