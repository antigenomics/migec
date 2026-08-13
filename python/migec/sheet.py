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

import re
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
                    # A positional pattern has nothing to score, so it must be matched the way it
                    # will be run -- anchored. Compiling it against a free scan raises the very
                    # error the anchor exists to avoid.
                    _core.match_pattern(
                        spec, "A" * (len(spec) + 1), max_offset=0 if is_positional(spec) else -1
                    )
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


# --------------------------------------------------------------------------------------------
# Positional layouts: the primary way to say where a barcode is.
#
# Most libraries put the barcode at a fixed offset in one read and have no constant sequence to
# anchor on. Two spellings, both anchored at the first base by construction:
#
#     ^NNNNNNNN          a pattern, pinned to position 0 by the leading caret
#     ^NNNNXXXXXXXX      N is a UMI base, X a cell barcode base
#     0:8                a half-open slice, the same 8 nt UMI
#     0:4,5:10           two slices: a 9 nt UMI split by one skipped base
#     cell:0:16,16:26    10x -- a 16 nt cell barcode then a 10 nt UMI
#
# Slices are half-open like Python's, which is the convention every downstream tool's docs use;
# `0:8` is eight bases and the next slice may start at 8. Gaps between slices become `.` (skipped,
# neither scored nor captured), which is exactly what a spacer is.

_SLICE = re.compile(r"^(?:(umi|cell):)?(\d+):(\d+)$")


def _from_slices(spec: str) -> str:
    """Translate `0:4,5:10` (or `cell:0:16,16:26`) into a pattern."""
    out: list[str] = []
    pos = 0
    for part in spec.split(","):
        m = _SLICE.match(part.strip())
        if not m:
            raise ValueError(
                f"{part.strip()!r} is not a slice -- expected START:STOP, half-open and 0-based, "
                f"optionally prefixed `umi:` or `cell:` (default umi). A 10 nt UMI at the read "
                f"start is `0:10`; 10x is `cell:0:16,16:26`."
            )
        kind, start, stop = m.group(1) or "umi", int(m.group(2)), int(m.group(3))
        if stop <= start:
            raise ValueError(
                f"{part.strip()!r}: slices are half-open, so STOP must exceed START -- "
                f"a 4 nt barcode at the read start is `0:4`, not `0:3`."
            )
        if start < pos:
            raise ValueError(
                f"{part.strip()!r} starts before the previous slice ended ({pos}). Slices must be "
                f"in increasing order and must not overlap: one base belongs to one barcode."
            )
        out.append("." * (start - pos))
        out.append(("N" if kind == "umi" else "X") * (stop - start))
        pos = stop
    return "".join(out)


def parse_layout(spec: str) -> tuple[str, bool]:
    """Translate a layout spec into `(pattern, anchored)`.

    `anchored` is True when the layout fixes the barcode at the first base, which is what
    `--max-offset 0` means. A caret says so explicitly; a slice list says so by construction,
    because a position is only a position if it is measured from somewhere.
    """
    text = spec.strip()
    anchored = text.startswith("^")
    if anchored:
        text = text[1:].strip()
    if not text:
        raise ValueError("empty layout")
    # A colon is not in the pattern grammar, so its presence is unambiguous.
    if ":" in text:
        return _from_slices(text), True
    return text, anchored


def is_positional(pattern: str) -> bool:
    """True when nothing in the pattern can be scored, so there is nothing to anchor a scan on.

    A free scan over such a pattern has no evidence to choose an offset with and `compile()`
    refuses it. That refusal is correct, so the offset is settled here instead of asking the user
    to pair every positional layout with `--max-offset 0` by hand.
    """
    return all(c in "NnXx." for c in pattern) and pattern != ""


# Layouts that come up often enough to be worth a name. Each is a real, published chemistry; the
# `source` is where the layout is written down, so a wrong one is falsifiable rather than folklore.
PRESETS: dict[str, tuple[str, str | None, str]] = {
    # name:          (master pattern,            slave pattern, description)
    "umi": (
        "^NNNNNNNN",
        None,
        "generic inline UMI, 8 nt at the start of the read. Change the run length to match yours "
        "or write the slice, `0:12`.",
    ),
    "migec": (
        "cagtggtatcaacgcagagtNNNNtNNNNtNNNN",
        None,
        "MIGEC 5'-RACE RepSeq: the SMART adapter then a 12 nt UMI split by two spacer bases. "
        "Source: misc/barcodes.txt of MIGEC 1.2.9 (tag v1-final). Prefix a sample tag such as "
        "`aaACT` per row to demultiplex.",
    ),
    "primerid": (
        "NNNNNNNNNcagtttaacttttgggccatcca",
        None,
        "HIV-1 Primer ID amplicon as used by MAGERI: a 9 nt UMI ahead of the gene-specific primer. "
        "Source: recovered by `migec suggest` from SRR1763769; the primer places it, so the scan "
        "stays free.",
    ),
    "duplex": (
        "^NNNNNNNNNNNN.....",
        "^NNNNNNNNNNNN.....",
        "Duplex sequencing (Schmitt/Kennedy): a 12 nt UMI and a 5 nt fixed spacer on BOTH mates. "
        "The two halves concatenate into one 24 nt strand-aware identifier. Warning: migec emits "
        "single-strand consensuses; pairing the two strands into a duplex consensus is not "
        "implemented, so do not quote a duplex error rate from this.",
    ),
    "10x": (
        "^XXXXXXXXXXXXXXXXNNNNNNNNNNNN",
        None,
        "10x Chromium 3' v3/v3.1: 16 nt cell barcode then a 12 nt UMI on R1. Run the later stages "
        "on R2 -- R1 is barcode and nothing else.",
    ),
    "10x-v2": (
        "^XXXXXXXXXXXXXXXXNNNNNNNNNN",
        None,
        "10x Chromium 3' v2 and 5' v1/v2: 16 nt cell barcode then a 10 nt UMI on R1.",
    ),
    "tso500": (
        "^NNNNN.....",
        None,
        "Illumina TSO500 ctDNA: a 5 nt UMI and a 5 nt spacer, on R1 ONLY -- the fgbio read "
        "structure is `5M5S+T +T`, R2 being all template. Warning: 5 nt is 1,024 barcodes, so the "
        "UMI alone does NOT identify a molecule on a real panel; TSO500's own pipeline groups on "
        "the mapping position as well (fgbio GroupReadsByUmi, after alignment). migec groups on "
        "the barcode, so it will report the space as saturated -- believe it.",
    ),
    "smarter-umi": (
        "^NNNNNNNNNNGGG",
        None,
        "SMARTer template-switching RNA-seq with a 10 nt inline UMI, then the GGG the template "
        "switch leaves behind. Source: ncgr/UMI-analysis, `fastq_qual_filter ... 0 10`.",
    ),
}


def preset(name: str) -> tuple[str, str | None]:
    """Look up a named layout. Raises ValueError listing every preset when the name is unknown."""
    key = name.strip().lower()
    if key not in PRESETS:
        listed = "\n".join(f"  {n:<12} {PRESETS[n][0]}" for n in PRESETS)
        raise ValueError(f"unknown preset {name!r}. Available:\n{listed}")
    master, slave, _ = PRESETS[key]
    return master, slave


def format_presets() -> str:
    """The preset table, for `migec sheet --presets` and the docs."""
    out = []
    for name, (master, slave, description) in PRESETS.items():
        out.append(f"{name}")
        out.append(f"    pattern  {master}" + (f"    slave  {slave}" if slave else ""))
        out.append(f"    {description}")
    return "\n".join(out)


# fgbio/Picard read-structure syntax, which is what TSO500, fgbio and samtools all speak.
_SEGMENT = re.compile(r"(\d+|\+)([MBTS])")


def from_read_structure(structure: str) -> str:
    """Translate an fgbio read structure into a migec pattern.

    `5M5S+T` is TSO500: five UMI bases, a five-base spacer, then template. `16B10M+T` is 10x.
    Both are positional, so both carry their own anchor -- the caller does not set an offset.

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
