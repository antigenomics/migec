"""The wider UMI search must preserve the released behavior of cell-tagged cohorts."""
from __future__ import annotations

import random

import pytest

pytest.importorskip("migec._core", reason="build the native extension first")
from migec.refine import run
from migec.checkout import run as checkout_run


def test_multi_substitution_extension_does_not_reinterpret_cell_tags(tmp_path):
    # Every tagged cell/UMI pair is a genuine molecule. Two added classes differ
    # from the abundant one in two cell bases or two UMI bases, with confident Q40.
    parent_cell = "ACGTACGTACGTACGT"
    parent_umi = "ACGTACGTACGTAC"
    groups = [(parent_cell, parent_umi, 20),
              ("TCGTATGTACGTACGT", parent_umi, 1),
              (parent_cell, "TCGTATGTACGTAC", 1)]
    seen = {(cell, umi) for cell, umi, _ in groups}
    rng = random.Random(145)
    while len(groups) < 2003:
        cell = "".join(rng.choices("ACGT", k=16))
        umi = "".join(rng.choices("ACGT", k=14))
        if (cell, umi) not in seen:
            seen.add((cell, umi))
            groups.append((cell, umi, 1))
    reads = tmp_path / "tagged.fq"
    with reads.open("w") as handle:
        n = 0
        for cell, umi, depth in groups:
            for _ in range(depth):
                handle.write(f"@r{n} CB:Z:{cell}\tRX:Z:{umi}\tCY:Z:{'I' * 16}\tQX:Z:{'I' * 14}\n"
                             f"{'ACGT' * 20}\n+\n{'I' * 80}\n")
                n += 1
    result = run(reads, tmp_path / "refine", sample_id="S1")
    assert result["cell_length"] == 16
    assert result["reads"] == 2022
    assert result["barcodes"] == 2003
    assert result["merged"] == 0
    assert result["molecules"] == 2003

    # Checkout's molecule diagnostics must use the same assay context and CPU budget.
    raw = tmp_path / "raw.fq"
    with raw.open("w") as handle:
        n = 0
        for cell, umi, depth in groups:
            for _ in range(depth):
                handle.write(f"@r{n}\n{cell}{umi}{'ACGT' * 20}\n+\n{'I' * 110}\n")
                n += 1
    sheet = tmp_path / "barcodes.tsv"
    sheet.write_text("S1\tcell:0:16,16:30\n")
    checkout = checkout_run(raw, sheet, tmp_path / "checkout", threads=1)
    sample, = checkout["samples"]
    assert checkout["threads"] == 1
    assert sample["cell_length"] == 16
    assert sample["molecules_observed"] == 2003
