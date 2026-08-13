"""The downstream contract: what a consensus record has to look like for anything to read it.

`docs/downstream.rst` records what each tool did when it was run against real output. This file is
the part that needs nothing installed -- the shape of the record itself -- plus the same aligner
and quantifier checks, skipped when the tool is not on PATH.

The contract has two halves and they are consumed separately: the read NAME carries the molecule
id for tools that drop FASTQ comments (arda goes through dnaio, which does), and the tab-separated
SAM tags carry everything for tools that keep it (`minimap2 -y`, `bwa mem -C`).
"""

from __future__ import annotations

import gzip
import random
import shutil
import subprocess

import pytest

from migec.assemble import run as assemble_run
from migec.checkout import run as checkout_run
from migec.refine import run as refine_run

CELLS, MOLECULES, READS = 6, 5, 3


@pytest.fixture(scope="module")
def consensus(tmp_path_factory):
    """A small 10x-shaped library, taken all the way to consensus FASTQ."""
    tmp = tmp_path_factory.mktemp("downstream")
    rng = random.Random(11)
    r1, r2 = tmp / "r1.fq.gz", tmp / "r2.fq.gz"
    with gzip.open(r1, "wt") as f1, gzip.open(r2, "wt") as f2:
        for c in range(CELLS):
            cb = "".join(rng.choice("ACGT") for _ in range(16))
            for m in range(MOLECULES):
                umi = "".join(rng.choice("ACGT") for _ in range(10))
                payload = "".join(rng.choice("ACGT") for _ in range(60))
                for r in range(READS):
                    f1.write(f"@c{c}m{m}r{r}\n{cb}{umi}\n+\n{'I' * 26}\n")
                    f2.write(f"@c{c}m{m}r{r}\n{payload}\n+\n{'I' * 60}\n")

    (tmp / "bc.txt").write_text("S\t" + "X" * 16 + "N" * 10 + "\n")
    checkout_run(r1, tmp / "bc.txt", tmp / "co", reads2=r2)
    refine_run(tmp / "co" / "S_R2.fq.gz", tmp / "ref", sample_id="S", expect_cells=CELLS)
    assemble_run(tmp / "ref" / "S.fq.gz", tmp / "asm", sample_id="S")
    return tmp / "asm" / "S.consensus.fq.gz"


def records(path):
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    return [lines[i : i + 4] for i in range(0, len(lines), 4)]


def test_one_record_is_one_molecule(consensus):
    assert len(records(consensus)) == CELLS * MOLECULES


def test_the_name_alone_identifies_the_molecule(consensus):
    """arda and STAR keep the name and drop the comment, so the name has to be self-sufficient."""
    names = set()
    for header, _, _, _ in records(consensus):
        name = header[1:].split(" ", 1)[0].split("\t", 1)[0]
        sample, cell, umi = name.split(".")
        assert (sample, len(cell), len(umi)) == ("S", 16, 10)
        names.add(name)
    assert len(names) == CELLS * MOLECULES  # and it is unique, or it is not an identifier


def test_the_tags_are_tab_separated_and_valid_sam(consensus):
    """Never: spaces here would produce a SAM line no parser accepts once an aligner appends it."""
    for header, _, _, _ in records(consensus):
        _, _, comment = header[1:].partition(" ")
        assert comment and "\t" in comment
        assert " " not in comment, "a space inside the comment splits one tag into two"
        tags = {}
        for field in comment.split("\t"):
            tag, kind, value = field.split(":", 2)
            assert len(tag) == 2 and kind in "ZifAB"
            tags[tag] = value
        assert set(tags) >= {"RX", "BC", "CB", "MI"}
        assert tags["MI"] == header[1:].split(" ", 1)[0]  # the name, again, for tools that tag


def _tags_in_sam(sam_text):
    counts: dict[str, int] = {}
    for line in sam_text.splitlines():
        if line.startswith("@"):
            continue
        for field in line.split("\t")[11:]:
            if field[:5] in ("RX:Z:", "CB:Z:", "MI:Z:"):
                counts[field[:2]] = counts.get(field[:2], 0) + 1
    return counts


@pytest.mark.parametrize("tool", ["minimap2", "bwa"])
def test_an_aligner_carries_the_tags_into_the_sam(consensus, tmp_path, tool):
    """The one property `-y` / `-C` exist for. Skipped where the aligner is not installed."""
    if shutil.which(tool) is None:
        pytest.skip(f"{tool} is not on PATH")

    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\n" + "".join(r[1] for r in records(consensus)) + "\n")
    n = len(records(consensus))

    if tool == "minimap2":
        cmd = ["minimap2", "-ax", "sr", "-y", str(reference), str(consensus)]
    else:
        subprocess.run(["bwa", "index", str(reference)], check=True, capture_output=True)
        cmd = ["bwa", "mem", "-C", str(reference), str(consensus)]
    sam = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout

    assert _tags_in_sam(sam) == {"RX": n, "CB": n, "MI": n}
