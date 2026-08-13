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
