"""Barcode metadata tables.

MIGEC's ``barcodes.txt`` is the format the published tables are written in, so it is read
verbatim -- tab or whitespace separated, ``#`` comments, one row per sample::

    S1<TAB>aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
    S2<TAB>aaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN

Column 3 onwards (slave barcode, R1/R2 paths) are accepted and ignored for now; paired-end
handling lands with the paired checkout. Several rows may share a sample id, which is how a
sample sequenced with more than one tag is declared.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SampleRow:
    sample_id: str
    pattern: str
    slave: str | None = None
    r1: str | None = None
    r2: str | None = None


def read_barcodes(path: str | Path) -> list[SampleRow]:
    """Parse a MIGEC-style barcode table. Raises ValueError with the line number on bad input."""
    rows: list[SampleRow] = []
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t") if "\t" in line else line.split()
            if len(fields) < 2:
                raise ValueError(
                    f"{path}:{lineno}: expected at least SAMPLE_ID and a pattern, got {line!r}"
                )
            sample_id, pattern = fields[0], fields[1]
            if not pattern.strip("."):
                raise ValueError(f"{path}:{lineno}: empty pattern for sample {sample_id!r}")

            def opt(i: int) -> str | None:
                if len(fields) <= i:
                    return None
                v = fields[i].strip()
                return None if v in ("", ".") else v

            rows.append(SampleRow(sample_id, pattern, opt(2), opt(3), opt(4)))

    if not rows:
        raise ValueError(f"{path}: no sample rows found")
    return rows


def describe(rows: list[SampleRow]) -> str:
    """A one-line-per-sample summary of what each pattern will extract."""
    from migec import _core

    out = []
    for r in rows:
        probe = _core.match_pattern(r.pattern, r.pattern.upper().replace("N", "A").replace(".", "A"))
        umi_len = len(probe["umi"])
        out.append(
            f"{r.sample_id}\tlen={len(r.pattern)}\tumi={umi_len}\tpattern={r.pattern}"
        )
    return "\n".join(out)
