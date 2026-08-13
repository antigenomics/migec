"""Barcode metadata tables.

MIGEC's ``barcodes.txt`` is the format the published tables are written in, so it is read
verbatim -- tab or whitespace separated, ``#`` comments, one row per sample::

    S1<TAB>aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
    S2<TAB>aaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN

Column 3 is MIGEC's *slave* barcode: a second pattern that sits on the other mate, whose captured
positions **extend** the UMI rather than starting a new one. That is how a 24 nt dual-end UMI is
declared::

    S1<TAB>NNNNNNNNNNNNtgact<TAB>agtcaNNNNNNNNNNNN

Both halves must match or the read is unmatched -- accepting the master alone would emit
half-length UMIs next to full-length ones, and every collision estimate downstream would then be
computed over two barcode spaces at once.

Columns 4-5 (R1/R2 paths) are accepted and ignored. Several rows may share a sample id, which is
how a sample sequenced with more than one tag is declared.
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
    """A one-line-per-sample summary of what each row will extract.

    Counts the captured positions from the spec rather than by matching a probe: `X` is not a
    IUPAC symbol, so a probe built by upper-casing the pattern carries an `X` base and the matcher
    has nothing sensible to do with it. Counting is also what the reader wants to check.
    """
    from migec import _core

    out = []
    for r in rows:
        # Compile every pattern, master and slave. `migec sheet` exists to catch a bad row before
        # a run rather than after one, so it has to actually parse them -- counting characters
        # would happily accept a symbol the grammar rejects.
        for spec in (r.pattern, r.slave):
            if spec:
                try:
                    _core.match_pattern(spec, "A" * (len(spec) + 1))
                except RuntimeError as exc:
                    raise ValueError(f"{r.sample_id}: {exc}") from exc
        umi = sum(1 for c in r.pattern if c in "Nn")
        cell = sum(1 for c in r.pattern if c in "Xx")
        # The slave pattern EXTENDS the UMI -- reporting only the master would tell a dual-end
        # sheet it has a 12 nt barcode when it has 24, which is the number every collision
        # estimate downstream is computed from.
        slave_umi = sum(1 for c in (r.slave or "") if c in "Nn")
        fields = [r.sample_id, f"len={len(r.pattern)}", f"umi={umi + slave_umi}"]
        if slave_umi:
            fields.append(f"({umi}+{slave_umi} dual-end)")
        if cell:
            fields.append(f"cell={cell}")
        fields.append(f"pattern={r.pattern}")
        if r.slave:
            fields.append(f"slave={r.slave}")
        out.append("\t".join(fields))
    return "\n".join(out)


# fgbio/Picard read-structure syntax, which is what TSO500, fgbio and samtools all speak.
_SEGMENT = __import__("re").compile(r"(\d+|\+)([MBTS])")


def from_read_structure(structure: str) -> str:
    """Translate an fgbio read structure into a migec pattern.

    `5M5S+T` is TSO500: five UMI bases, a five-base spacer, then template. `16B10M+T` is 10x.
    Both are positional, so both need `--max-offset 0`.

        M  molecular barcode   -> N   captured as UMI
        B  sample barcode      -> X   captured as cell barcode
        S  skip                -> .   neither scored nor captured
        T  template            -> the payload; the pattern ends here

    The pattern stops at the first template segment, because everything after it is the read and
    migec trims to exactly that point. `+` means "the rest of the read" and is only meaningful on
    the last segment.
    """
    text = structure.strip().upper()
    if not text:
        raise ValueError("empty read structure")
    segments = _SEGMENT.findall(text)
    if not segments or "".join(a + b for a, b in segments) != text:
        raise ValueError(
            f"{structure!r} is not an fgbio read structure -- expected runs like 5M5S+T, "
            f"where M is a UMI base, B a sample/cell barcode, S a skipped base and T template"
        )
    out = []
    for i, (count, kind) in enumerate(segments):
        if count == "+":
            if i != len(segments) - 1:
                raise ValueError(f"{structure!r}: '+' is only valid on the last segment")
            if kind != "T":
                raise ValueError(
                    f"{structure!r}: '+{kind}' captures an unbounded barcode, which has no length "
                    f"for the collision arithmetic to use. Give it a count."
                )
            break
        if kind == "T":
            break
        out.append({"M": "N", "B": "X", "S": "."}[kind] * int(count))
    pattern = "".join(out)
    if not pattern:
        raise ValueError(f"{structure!r} captures nothing before the template")
    return pattern
