"""The barcode table parser and the command line.

Both are trust boundaries with nothing behind them: a bad sheet is user input, and a CLI command
that does not run is a typo nobody catches until someone types it. Neither had a test.
"""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

from migec.cli import app
from migec.sheet import describe, read_barcodes

runner = CliRunner()


def sheet(tmp_path, text):
    p = tmp_path / "bc.txt"
    p.write_text(text)
    return p


# ------------------------------------------------------------------ the parser


def test_a_plain_table_parses(tmp_path):
    rows = read_barcodes(sheet(tmp_path, "S1\taaACTcagtggNNNNtNNNNtNNNN\nS2\taaAGAcagtggNNNNtNNNNtNNNN\n"))
    assert [r.sample_id for r in rows] == ["S1", "S2"]
    assert rows[0].slave is None


def test_whitespace_and_comments(tmp_path):
    rows = read_barcodes(sheet(tmp_path, "# a comment\n\nS1   ACGTNNNN\n"))
    assert len(rows) == 1
    assert rows[0].pattern == "ACGTNNNN"


def test_the_slave_column_is_read(tmp_path):
    rows = read_barcodes(sheet(tmp_path, "S1\tNNNNNNNNNNNNTGACT\tAGTCANNNNNNNNNNNN\n"))
    assert rows[0].slave == "AGTCANNNNNNNNNNNN"
    # A "." means absent, which is how a sheet declares a slave for some rows and not others.
    rows = read_barcodes(sheet(tmp_path, "S1\tACGTNNNN\t.\n"))
    assert rows[0].slave is None


def test_bad_input_names_the_line(tmp_path):
    with pytest.raises(ValueError, match=r":2:"):
        read_barcodes(sheet(tmp_path, "S1\tACGTNNNN\nS2\n"))
    with pytest.raises(ValueError, match="empty pattern"):
        read_barcodes(sheet(tmp_path, "S1\t...\n"))
    with pytest.raises(ValueError, match="no sample rows"):
        read_barcodes(sheet(tmp_path, "# nothing but a comment\n"))


def test_describe_counts_every_captured_position(tmp_path):
    """`migec sheet` answers "what will this row extract?", so it has to count the whole barcode."""
    rows = read_barcodes(
        sheet(tmp_path,
              "S1\taaACTcagNNNNtNNNNtNNNN\n"
              "S2\tNNNNNNNNNNNNTGACT\tAGTCANNNNNNNNNNNN\n"
              "S3\tXXXXXXXXXXXXXXXXNNNNNNNNNNNNcagtgg\n")
    )
    lines = describe(rows).splitlines()
    assert "umi=12" in lines[0]
    # Dual-end: both halves, and the split shown so the sheet can be checked against the protocol.
    assert "umi=24" in lines[1] and "12+12 dual-end" in lines[1]
    assert "slave=AGTCANNNNNNNNNNNN" in lines[1]
    # Cell barcode positions are captured separately and must not be lost in the UMI count.
    assert "umi=12" in lines[2] and "cell=16" in lines[2]


# ------------------------------------------------------------------ the command line


def test_every_command_runs(tmp_path):
    """Seven commands, none of which had a test that it starts."""
    assert runner.invoke(app, ["--help"]).exit_code == 0
    assert runner.invoke(app, ["info"]).exit_code == 0
    for cmd in ("checkout", "suggest", "refine", "assemble", "subsample", "sheet", "info"):
        result = runner.invoke(app, [cmd, "--help"])
        assert result.exit_code == 0, f"`migec {cmd} --help` failed:\n{result.output}"


def test_sheet_command_prints_the_table(tmp_path):
    p = sheet(tmp_path, "S1\taaACTcagNNNNtNNNNtNNNN\n")
    result = runner.invoke(app, ["sheet", str(p)])
    assert result.exit_code == 0
    assert "umi=12" in result.output


def test_a_bad_trim_mode_is_rejected(tmp_path):
    p = sheet(tmp_path, "S1\tACGTNNNN\n")
    result = runner.invoke(
        app, ["checkout", "reads.fq", "-b", str(p), "-o", str(tmp_path), "--trim", "sideways"]
    )
    assert result.exit_code != 0


def test_the_pipeline_runs_end_to_end_through_the_cli(tmp_path):
    """checkout -> refine -> assemble, driven the way a user drives it."""
    import gzip
    import random

    rng = random.Random(0)
    adapter = "CAGTGGTATCAACGCAGAGT"
    with gzip.open(tmp_path / "r.fq.gz", "wt") as fh:
        for i in range(400):
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            payload = "".join(rng.choice("ACGT") for _ in range(60))
            for _ in range(4):
                s = umi + adapter + payload
                fh.write(f"@r{i}_{_}\n{s}\n+\n{'I' * len(s)}\n")
    p = sheet(tmp_path, f"S1\t{'N' * 12}{adapter.lower()}\n")

    r = runner.invoke(app, ["checkout", str(tmp_path / "r.fq.gz"), "-b", str(p),
                            "-o", str(tmp_path / "co")])
    assert r.exit_code == 0, r.output
    r = runner.invoke(app, ["refine", str(tmp_path / "co" / "S1.fq.gz"),
                            "-o", str(tmp_path / "ref")])
    assert r.exit_code == 0, r.output
    r = runner.invoke(app, ["assemble", str(tmp_path / "ref" / "S1.fq.gz"),
                            "-o", str(tmp_path / "asm")])
    assert r.exit_code == 0, r.output
    assert (tmp_path / "asm" / "S1.consensus.fq.gz").exists()

    r = runner.invoke(app, ["subsample", str(tmp_path / "co" / "S1.fq.gz"),
                            "-o", str(tmp_path / "sub.fq.gz"), "--keep", "25"])
    assert r.exit_code == 0, r.output
    assert (tmp_path / "sub.fq.gz").exists()


def test_describe_rejects_a_pattern_the_grammar_does_not_accept(tmp_path):
    """`migec sheet` exists to catch a bad row before a run, not after one."""
    rows = read_barcodes(sheet(tmp_path, "S1\tACGT?NNN\n"))
    with pytest.raises(ValueError, match="S1:"):
        describe(rows)
    # ...including in the slave column, which is just as capable of being wrong.
    rows = read_barcodes(sheet(tmp_path, "S1\tACGTNNNN\tAGTC?NNN\n"))
    with pytest.raises(ValueError, match="S1:"):
        describe(rows)
    # A pattern with no scored position matches everywhere and is refused by the compiler.
    rows = read_barcodes(sheet(tmp_path, "S1\tNNNN\n"))
    with pytest.raises(ValueError, match="S1:"):
        describe(rows)
