"""The barcode table parser and the command line.

Both are trust boundaries with nothing behind them: a bad sheet is user input, and a CLI command
that does not run is a typo nobody catches until someone types it. Neither had a test.
"""

from __future__ import annotations

import re

import pytest
from typer.testing import CliRunner

from migec.cli import app
from migec.sheet import describe, read_barcodes

# A wide terminal, because typer wraps its error box to the terminal width and CI runs at 80
# columns: an assertion on a message can otherwise pass locally and fail there, having found the
# message split across two lines. `clean()` is the belt to that braces -- it also strips the ANSI
# colour codes typer emits.
runner = CliRunner(env={"COLUMNS": "200", "TERM": "dumb"})

_ANSI = re.compile(r"\x1b\[[0-9;]*m")


def clean(output: str) -> str:
    """Output with colour codes gone and whitespace collapsed, so a wrapped line still matches."""
    return " ".join(_ANSI.sub("", output).split())


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
    assert "umi=12" in clean(result.output)


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
    # A pattern with no scored position has nowhere to be scanned for -- but it does not need to
    # be, because it is positional, so `describe` matches it the way a run would: anchored.
    # Refusing it here is what used to force `--max-offset 0` into every 10x command line.
    rows = read_barcodes(sheet(tmp_path, "S1\tNNNN\n"))
    assert "umi=4" in describe(rows)


# ------------------------------------------------------- inline patterns (umi_tools style)


def test_bc_pattern_replaces_the_sheet(tmp_path):
    """A positional chemistry is one sample and one pattern; writing a one-line file for it is
    friction umi_tools, umitools and mgatk all avoid."""
    import gzip

    with gzip.open(tmp_path / "r1.fq.gz", "wt") as f1, gzip.open(tmp_path / "r2.fq.gz", "wt") as f2:
        for i in range(50):
            f1.write(f"@r{i}\n{'ACGT' * 4}{'TTTTTTTTTT'}\n+\n{'I' * 26}\n")
            f2.write(f"@r{i}\n{'G' * 90}\n+\n{'I' * 90}\n")
    r = runner.invoke(app, ["checkout", str(tmp_path / "r1.fq.gz"), str(tmp_path / "r2.fq.gz"),
                            "--bc-pattern", "X" * 16 + "N" * 10, "--sample", "PBMC",
                            "--max-offset", "0", "-o", str(tmp_path / "out")])
    assert r.exit_code == 0, r.output
    assert (tmp_path / "out" / "PBMC_R1.fq.gz").exists()


def test_a_umi_tools_pattern_is_refused_rather_than_misread(tmp_path):
    """Never: umi_tools spells a cell barcode `C`; here `C` is cytosine. Pasting one would compile
    into a pattern demanding literal cytosines, match nothing, and look like a bad library."""
    r = runner.invoke(app, ["checkout", "reads.fq", "--bc-pattern", "CCCCCCCCCCCCCCCCNNNNNNNNNN",
                            "-o", str(tmp_path / "out")])
    assert r.exit_code != 0
    assert "cytosine" in clean(r.output)
    assert "XXXXXXXXXXXXXXXXNNNNNNNNNN" in clean(r.output)


def test_exactly_one_of_barcodes_and_bc_pattern(tmp_path):
    p = sheet(tmp_path, "S1\tACGTNNNN\n")
    for args in ([], ["-b", str(p), "--bc-pattern", "ACGTNNNN"]):
        r = runner.invoke(app, ["checkout", "reads.fq", "-o", str(tmp_path / "o"), *args])
        assert r.exit_code != 0
        assert "exactly one" in clean(r.output)


# ------------------------------------------------- fgbio read structures (TSO500, 10x)


def test_read_structure_translates_the_platforms(tmp_path):
    from migec.sheet import from_read_structure

    assert from_read_structure("5M5S+T") == "NNNNN....."          # TSO500
    assert from_read_structure("16B10M+T") == "X" * 16 + "N" * 10  # 10x 5'
    assert from_read_structure("8M+T") == "N" * 8
    assert from_read_structure("5m5s+t") == "NNNNN....."           # case-insensitive


def test_read_structure_refuses_what_it_cannot_size(tmp_path):
    from migec.sheet import from_read_structure

    # An unbounded barcode has no length for the collision arithmetic to use.
    with pytest.raises(ValueError, match="unbounded"):
        from_read_structure("5M+M")
    with pytest.raises(ValueError, match="last segment"):
        from_read_structure("+T5M")
    with pytest.raises(ValueError, match="fgbio read structure"):
        from_read_structure("5Z")
    with pytest.raises(ValueError, match="captures nothing"):
        from_read_structure("+T")


def test_tso500_read_structure_end_to_end(tmp_path):
    """TSO500: 5 nt UMI, 5 nt spacer, then template -- on both mates, so the two halves
    concatenate into one 10 nt molecule identifier."""
    import gzip
    import random

    rng = random.Random(0)
    truth = []
    with gzip.open(tmp_path / "r1.fq.gz", "wt") as f1, gzip.open(tmp_path / "r2.fq.gz", "wt") as f2:
        for i in range(300):
            u1 = "".join(rng.choice("ACGT") for _ in range(5))
            u2 = "".join(rng.choice("ACGT") for _ in range(5))
            truth.append(u1 + u2)
            s1 = u1 + "ACGTA" + "".join(rng.choice("ACGT") for _ in range(80))
            s2 = u2 + "ACGTA" + "".join(rng.choice("ACGT") for _ in range(80))
            f1.write(f"@r{i}\n{s1}\n+\n{'I' * len(s1)}\n")
            f2.write(f"@r{i}\n{s2}\n+\n{'I' * len(s2)}\n")

    r = runner.invoke(app, ["checkout", str(tmp_path / "r1.fq.gz"), str(tmp_path / "r2.fq.gz"),
                            "--read-structure", "5M5S+T", "--read-structure2", "5M5S+T",
                            "--sample", "TSO", "--max-offset", "0", "-o", str(tmp_path / "out")])
    assert r.exit_code == 0, r.output

    seen = []
    with gzip.open(tmp_path / "out" / "TSO_R1.fq.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                seen.append(next(f[5:] for f in line.split() if f.startswith("RX:Z:")))
    assert len(seen) == 300
    assert all(len(u) == 10 for u in seen), "both mates' UMIs, concatenated"
    assert set(seen) == set(truth)


def test_read_structure_and_bc_pattern_are_exclusive(tmp_path):
    r = runner.invoke(app, ["checkout", "reads.fq", "--read-structure", "5M5S+T",
                            "--bc-pattern", "NNNNN", "-o", str(tmp_path / "o")])
    assert r.exit_code != 0 and "exactly one of" in clean(r.output)
    r = runner.invoke(app, ["checkout", "reads.fq", "--read-structure2", "5M5S+T",
                            "-o", str(tmp_path / "o")])
    assert r.exit_code != 0 and "needs --read-structure" in clean(r.output)


# ------------------------------------------------------------------ positional layouts

def test_the_two_positional_spellings_agree():
    """`^NNNN...` and a slice list must compile to the same pattern, or the docs lie."""
    from migec.sheet import parse_layout

    assert parse_layout("^NNNNNNNN") == ("NNNNNNNN", True)
    assert parse_layout("0:8") == ("NNNNNNNN", True)
    assert parse_layout("0:4,5:10") == ("NNNN.NNNNN", True)
    assert parse_layout("cell:0:16,16:26") == ("X" * 16 + "N" * 10, True)


def test_a_pattern_with_an_anchor_is_not_anchored_by_default():
    """An adapter places the pattern, so the scan stays free -- that is the whole point of it."""
    from migec.sheet import parse_layout

    assert parse_layout("NNNNcagtggtatcaacgcagagt") == ("NNNNcagtggtatcaacgcagagt", False)


def test_half_open_slices_are_rejected_when_empty_or_out_of_order():
    from migec.sheet import parse_layout

    for bad in ("0:0", "5:2"):
        with pytest.raises(ValueError, match="half-open"):
            parse_layout(bad)
    with pytest.raises(ValueError, match="before the previous slice"):
        parse_layout("0:8,4:12")
    with pytest.raises(ValueError, match="not a slice"):
        parse_layout("0:4,junk")


def test_is_positional_only_when_nothing_can_be_scored():
    from migec.sheet import is_positional

    assert is_positional("NNNNXXXX...")
    assert not is_positional("NNNNcagt")
    assert not is_positional("NNNNNGGG")


def test_every_preset_compiles_and_is_described():
    """A preset nobody can run is worse than no preset: it looks supported."""
    from migec.sheet import PRESETS, SampleRow, describe, parse_layout, preset

    for name in PRESETS:
        master, slave = preset(name)
        pattern, _ = parse_layout(master)
        rows = [SampleRow(name, pattern, parse_layout(slave)[0] if slave else None)]
        assert describe(rows)          # compiles both halves, or raises
        assert "N" in pattern          # every layout captures a UMI


def test_an_unknown_preset_lists_the_real_ones():
    from migec.sheet import PRESETS, preset

    with pytest.raises(ValueError) as exc:
        preset("nope")
    for name in PRESETS:
        assert name in str(exc.value)


def test_positional_needs_no_max_offset_on_the_command_line(tmp_path):
    """The regression this exists for: `--preset 10x-v2` used to need `--max-offset 0` spelled out,
    and without it checkout refused every read with an error about anchoring."""
    import gzip
    import random

    rng = random.Random(3)
    r1, r2 = tmp_path / "r1.fq.gz", tmp_path / "r2.fq.gz"
    with gzip.open(r1, "wt") as f1, gzip.open(r2, "wt") as f2:
        for c in range(5):
            cb = "".join(rng.choice("ACGT") for _ in range(16))
            for m in range(4):
                umi = "".join(rng.choice("ACGT") for _ in range(10))
                pay = "".join(rng.choice("ACGT") for _ in range(50))
                for r in range(3):
                    f1.write(f"@c{c}m{m}r{r}\n{cb}{umi}\n+\n{'I' * 26}\n")
                    f2.write(f"@c{c}m{m}r{r}\n{pay}\n+\n{'I' * 50}\n")

    for spec in (["--preset", "10x-v2"],
                 ["--bc-pattern", "cell:0:16,16:26"],
                 ["--bc-pattern", "^" + "X" * 16 + "N" * 10],
                 ["--read-structure", "16B10M+T"]):
        out = tmp_path / ("co" + spec[-1][:6].replace(":", "_"))
        result = runner.invoke(app, ["checkout", str(r1), str(r2), *spec, "-o", str(out)])
        assert result.exit_code == 0, clean(result.output)
        # 5 cells x 4 molecules, all of them found -- the number is the point, not that it ran.
        assert "20" in clean(result.output)


def test_giving_two_layouts_at_once_is_refused(tmp_path):
    result = runner.invoke(
        app, ["checkout", "x.fq", "--preset", "10x", "--bc-pattern", "0:8", "-o", str(tmp_path)]
    )
    assert result.exit_code != 0
    assert "exactly one of" in clean(result.output)


def test_sheet_presets_lists_them():
    from migec.sheet import PRESETS

    result = runner.invoke(app, ["sheet", "--presets"])
    assert result.exit_code == 0
    for name in PRESETS:
        assert name in result.output
