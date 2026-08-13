"""Dual-end barcodes: half the UMI on each mate.

MAGERI's design, quoted verbatim from its Methods: ``NNNNNNNNNNNNtgact`` on one mate and
``agtcaNNNNNNNNNNNN`` on the other, giving a 24 nt UMI from twelve bases at each end of the
molecule.
"""

from __future__ import annotations

import gzip
import random

import pytest

from migec.checkout import run

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
