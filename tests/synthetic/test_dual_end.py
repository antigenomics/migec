"""Dual-end barcodes: half the UMI on each mate.

MAGERI's design, quoted verbatim from its Methods: ``NNNNNNNNNNNNtgact`` on one mate and
``agtcaNNNNNNNNNNNN`` on the other, giving a 24 nt UMI from twelve bases at each end of the
molecule.
"""

from __future__ import annotations

import gzip
import random

import pytest

# Never: BEFORE the `migec` imports below, and not `pytestmark` alone. A module-scope
# import of a package whose extension is missing raises at COLLECTION, which pytest
# reports as an error rather than a skip -- so a machine without the built extension
# fails the suite instead of saying it cannot run it.
pytest.importorskip("migec._core", reason="the C++ extension is not built: run `bash setup.sh`")

from migec.checkout import run

from tests.conftest import requires_core

pytestmark = requires_core

MASTER = "NNNNNNNNNNNNTGACT"
SLAVE = "AGTCANNNNNNNNNNNN"


def corpus(tmp_path, n=2000, seed=0, break_slave=False):
    rng = random.Random(seed)
    truth = []
    with gzip.open(tmp_path / "r1.fq.gz", "wt") as f1, gzip.open(tmp_path / "r2.fq.gz", "wt") as f2:
        for i in range(n):
            u1 = "".join(rng.choice("ACGT") for _ in range(12))
            u2 = "".join(rng.choice("ACGT") for _ in range(12))
            truth.append(u1 + u2)
            s1 = u1 + "TGACT" + "".join(rng.choice("ACGT") for _ in range(60))
            handle = "TTTTT" if break_slave and i % 2 == 0 else "AGTCA"
            s2 = handle + u2 + "".join(rng.choice("ACGT") for _ in range(60))
            f1.write(f"@r{i}\n{s1}\n+\n{'I' * len(s1)}\n")
            f2.write(f"@r{i}\n{s2}\n+\n{'I' * len(s2)}\n")
    return truth


def test_the_umi_is_both_halves(tmp_path):
    truth = corpus(tmp_path)
    (tmp_path / "bc.txt").write_text(f"S1\t{MASTER}\t{SLAVE}\n")
    s = run(tmp_path / "r1.fq.gz", tmp_path / "bc.txt", tmp_path / "out",
            reads2=tmp_path / "r2.fq.gz", max_offset=0)
    assert s["assigned"] == 2000
    assert s["samples"][0]["umi_length"] == 24

    seen = []
    with gzip.open(tmp_path / "out" / "S1_R1.fq.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                seen.append(next(f[5:] for f in line.split() if f.startswith("RX:Z:")))
    assert len(seen) == 2000
    assert all(len(u) == 24 for u in seen)
    assert set(seen) == set(truth)


def test_both_halves_or_nothing(tmp_path):
    """Accepting the master alone would emit 12 nt UMIs beside 24 nt ones, and every collision
    estimate downstream would then be computed over two barcode spaces at once."""
    corpus(tmp_path, break_slave=True)
    (tmp_path / "bc.txt").write_text(f"S1\t{MASTER}\t{SLAVE}\n")
    s = run(tmp_path / "r1.fq.gz", tmp_path / "bc.txt", tmp_path / "out",
            reads2=tmp_path / "r2.fq.gz", max_offset=0)
    assert s["assigned"] == 1000
    assert s["unmatched"] == 1000

    with gzip.open(tmp_path / "out" / "S1_R1.fq.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                assert len(next(f[5:] for f in line.split() if f.startswith("RX:Z:"))) == 24


def test_a_slave_without_a_second_mate_is_refused(tmp_path):
    """A dual-end sheet run single-end cannot silently produce half a barcode."""
    corpus(tmp_path, n=10)
    (tmp_path / "bc.txt").write_text(f"S1\t{MASTER}\t{SLAVE}\n")
    with pytest.raises(RuntimeError, match="second mate"):
        run(tmp_path / "r1.fq.gz", tmp_path / "bc.txt", tmp_path / "out", max_offset=0)


def test_a_free_scan_correctly_refuses_a_five_base_handle(tmp_path):
    """TGACT occurs by chance about every kilobase, so an unanchored scan cannot place it -- and
    saying so is right. The fix is --max-offset, not a lower bar."""
    corpus(tmp_path, n=200)
    (tmp_path / "bc.txt").write_text(f"S1\t{MASTER}\t{SLAVE}\n")
    s = run(tmp_path / "r1.fq.gz", tmp_path / "bc.txt", tmp_path / "out",
            reads2=tmp_path / "r2.fq.gz")
    assert s["assigned"] == 0


# ------------------------------------------------------- positional chemistries (10x)

TENX = "X" * 16 + "N" * 10


def test_a_positional_barcode_read_with_no_payload_is_not_dropped(tmp_path):
    """10x: R1 is 26 nt of cell barcode and UMI and nothing else, R2 is the cDNA. Checking R1's
    leftover length by itself drops 100% of a perfectly good library as 'too short'."""
    with gzip.open(tmp_path / "r1.fq.gz", "wt") as f1, gzip.open(tmp_path / "r2.fq.gz", "wt") as f2:
        rng = random.Random(0)
        for i in range(500):
            cb = "".join(rng.choice("ACGT") for _ in range(16))
            umi = "".join(rng.choice("ACGT") for _ in range(10))
            cdna = "".join(rng.choice("ACGT") for _ in range(90))
            f1.write(f"@r{i}\n{cb}{umi}\n+\n{'I' * 26}\n")
            f2.write(f"@r{i}\n{cdna}\n+\n{'I' * 90}\n")
    (tmp_path / "bc.txt").write_text(f"P\t{TENX}\n")

    s = run(tmp_path / "r1.fq.gz", tmp_path / "bc.txt", tmp_path / "out",
            reads2=tmp_path / "r2.fq.gz", max_offset=0)
    assert s["assigned"] == 500
    assert s["short_payload"] == 0
    assert s["samples"][0]["umi_length"] == 10

    # Both mates carry the tags, and the cell barcode is one of them.
    with gzip.open(tmp_path / "out" / "P_R2.fq.gz", "rt") as fh:
        header = fh.readline()
    tags = dict(f.split(":Z:") for f in header.split() if ":Z:" in f)
    assert len(tags["CB"]) == 16
    assert len(tags["RX"]) == 10


def test_the_skew_warning_is_weighed_against_the_barcode_not_the_umi(tmp_path):
    """The counter is keyed on cell + UMI, so `effective_length` spans 26 nt, not 10.

    Weighed against `umi_length` the test reads "25.4 nt of a 10 nt UMI is usable" and can never
    be true, so the warning is dead on every single-cell run -- which is exactly where a skewed
    cell barcode is worth hearing about. Here the first eight positions of the cell barcode are
    constant, so about 18 of 26 barcode bases carry information: below 0.8 x 26, above 0.8 x 10.
    """
    with gzip.open(tmp_path / "r1.fq.gz", "wt") as f1, gzip.open(tmp_path / "r2.fq.gz", "wt") as f2:
        rng = random.Random(3)
        for i in range(500):
            cb = "A" * 8 + "".join(rng.choice("ACGT") for _ in range(8))
            umi = "".join(rng.choice("ACGT") for _ in range(10))
            cdna = "".join(rng.choice("ACGT") for _ in range(90))
            f1.write(f"@r{i}\n{cb}{umi}\n+\n{'I' * 26}\n")
            f2.write(f"@r{i}\n{cdna}\n+\n{'I' * 90}\n")
    (tmp_path / "bc.txt").write_text(f"P\t{TENX}\n")

    s = run(tmp_path / "r1.fq.gz", tmp_path / "bc.txt", tmp_path / "out",
            reads2=tmp_path / "r2.fq.gz", max_offset=0)["samples"][0]
    # The three lengths are reported apart, so nothing has to guess which one a column means.
    assert (s["umi_length"], s["cell_length"], s["barcode_length"]) == (10, 16, 26)
    assert s["umi_length"] < s["effective_length"] < s["barcode_length"]

    from migec.checkout import format_report

    report = format_report(run(tmp_path / "r1.fq.gz", tmp_path / "bc.txt", tmp_path / "out2",
                               reads2=tmp_path / "r2.fq.gz", max_offset=0))
    assert "barcode is 26 nt but only" in report


def test_a_cell_that_fills_the_key_on_its_own_does_not_shift_by_64(tmp_path):
    """32 nt of cell and no UMI: the one layout that reaches a 64-bit shift of a 64-bit key.

    `>> 64` is undefined, not zero -- the count is masked to 6 bits on x86, so it becomes `>> 0`
    and the UMI is ORed back on top of the cell barcode instead of below it. Verified under UBSan
    that the unguarded expression is flagged and returns a different value from the guarded one.
    The length check permits this layout (32 + 0 fits the key exactly), so the guard is what makes
    it defined: with no UMI the molecule key is the cell key, which is what these ten reads say.
    """
    with gzip.open(tmp_path / "r1.fq.gz", "wt") as f1, gzip.open(tmp_path / "r2.fq.gz", "wt") as f2:
        for i in range(10):
            f1.write(f"@r{i}\n{'ACGT' * 8}\n+\n{'I' * 32}\n")
            f2.write(f"@r{i}\n{'A' * 90}\n+\n{'I' * 90}\n")
    (tmp_path / "bc.txt").write_text(f"P\t{'X' * 32}\n")

    s = run(tmp_path / "r1.fq.gz", tmp_path / "bc.txt", tmp_path / "out",
            reads2=tmp_path / "r2.fq.gz", max_offset=0)["samples"][0]
    assert (s["cell_length"], s["umi_length"], s["barcode_length"]) == (32, 0, 32)
    assert s["reads"] == 10
    assert s["umis"] == 1, "ten reads of one cell barcode are one key, not ten"


def test_a_barcode_too_long_for_the_key_is_refused_by_name(tmp_path):
    """20 nt of cell and 16 of UMI clear their own bounds and their 36 nt sum does not fit at all.

    Before the guard this counted 1 distinct barcode where there were 2: the key holds 32 bases,
    so two molecules sharing their first 32 and differing past them packed identically and merged.
    That is the error nothing downstream can detect, which is why it is refused rather than warned
    about -- and refused when the pattern compiles, where the offending row can still be named.
    """
    (tmp_path / "bc.txt").write_text(f"P\t{'X' * 20}{'N' * 16}\n")
    with gzip.open(tmp_path / "r1.fq.gz", "wt") as f1, gzip.open(tmp_path / "r2.fq.gz", "wt") as f2:
        f1.write(f"@r0\n{'ACGT' * 9}\n+\n{'I' * 36}\n")
        f2.write(f"@r0\n{'A' * 90}\n+\n{'I' * 90}\n")
    with pytest.raises(RuntimeError, match="32 bases in total"):
        run(tmp_path / "r1.fq.gz", tmp_path / "bc.txt", tmp_path / "out",
            reads2=tmp_path / "r2.fq.gz", max_offset=0)


def test_a_positional_pattern_is_refused_by_a_free_scan(tmp_path):
    """No anchor means nothing to search for. Saying so beats matching everywhere.

    The refusal is still there, but you now have to ask for the free scan: a pattern with nothing
    to score is anchored by default, which is what the caret and the slice syntax say out loud.
    """
    with gzip.open(tmp_path / "r1.fq.gz", "wt") as f1, gzip.open(tmp_path / "r2.fq.gz", "wt") as f2:
        f1.write(f"@r0\n{'A' * 26}\n+\n{'I' * 26}\n")
        f2.write(f"@r0\n{'C' * 90}\n+\n{'I' * 90}\n")
    (tmp_path / "bc.txt").write_text(f"P\t{TENX}\n")
    with pytest.raises(RuntimeError, match="max_offset"):
        run(tmp_path / "r1.fq.gz", tmp_path / "bc.txt", tmp_path / "out",
            reads2=tmp_path / "r2.fq.gz", max_offset=-1)

    # ...and with the default it simply works, which is the point of the change.
    s = run(tmp_path / "r1.fq.gz", tmp_path / "bc.txt", tmp_path / "out2",
            reads2=tmp_path / "r2.fq.gz")
    assert s["assigned"] == 1
