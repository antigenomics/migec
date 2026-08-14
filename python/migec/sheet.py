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
        "^NNNNNNNNNN...",
        None,
        "SMARTer template-switching RNA-seq with a 10 nt inline UMI, then the GGG the template "
        "switch leaves behind. Source: ncgr/UMI-analysis, `fastq_qual_filter ... 0 10` and "
        "`fastq_umi_clipper <fq> 10 3`. Note: the GGG is skipped, not scored -- three bases are "
        "6.0 bits against an anchored bar of 6.64, so scoring them refuses every read. The source "
        "pipeline erases them without checking them either.",
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


# --------------------------------------------------------------------------------------------
# Assay profiles: what to run, once you know what the experiment is FOR.
#
# A layout says where the barcode is. It does not say how many reads a consensus needs before it is
# worth anything, and that turns out to matter more than the layout does. The two axes are
# independent -- the same 12 nt inline UMI serves a repertoire census and an MRD assay, and the
# right settings are opposite.
#
# Never: `--min-reads` defaults to 1, which is right for COUNTING molecules and wrong for CALLING
# VARIANTS. A consensus over one read IS that read -- no error correction at all, just counting.
# Measured on certified cfDNA reference material at 20 ng/10x, every systematic `-> G` false
# positive (the 2-colour dark-G artifact) vanished at `--min-reads 3`, taking specificity from 0%
# to 100%, while every certified variant survived with its frequency stable to the third decimal.
# It cost 41% of the molecules and bought an 8x better limit. `docs/detection.rst`.
#
# Never: a counting assay must NOT inherit the variant-calling threshold. Discarding singletons
# throws away 79% of the barcodes on a shallow repertoire library, and a clonotype seen once is
# still a clonotype. The `sensitivity` field is what separates them, and it is the whole point of
# this table:
#
#     counting        every molecule counts, errors are corrected not excluded   --min-reads 1
#     sensitive       a variant needs a real family behind it                    --min-reads 2
#     ultrasensitive  a variant needs a family that outvoted the chemistry       --min-reads 3
#
# Note: `payload_diverse` records whether the template itself distinguishes molecules. On an AIRR
# amplicon or an MRD IGH assay the rearrangement is near-unique, so `refine`'s payload-agreement
# term and `assemble`'s linkage sub-clustering both have real evidence to work with. On an exome,
# an RNA-seq library or a ctDNA panel every molecule at a locus reads the same, so that evidence is
# worth ~nothing and the barcode carries the whole burden. Never: payload agreement must be
# discounted by the measured clonality -- on a clonal library it is worth nothing even when the
# payload is nominally diverse.


# Never: the pre-amplification floor is NOT a reverse-transcription rate outside an RNA assay, and
# three of the seven profiles below are DNA. What the floor actually is, per class:
#
#   RNA (airr, rnaseq, 10x)   reverse transcriptase miscall, then the first PCR cycle. The RT is
#                             the dominant term and 1e-4 is 10x's own figure for their V(D)J RT.
#   DNA (exome, amplicon,     library-prep damage plus the first PCR cycle. Oxidation of guanine
#        ctdna, mrd)          during acoustic shearing gives 8-oxoG and therefore G>T / C>A;
#                             cytosine deamination gives C>T / G>A. Both are in the molecule
#                             before amplification, so every read of the group carries them and no
#                             consensus removes them -- which is the same argument as for RT, with
#                             a different chemistry supplying it.
#
# Note: that damage signature is NOT the artifact measured here. Ours is `-> G`, which is the
# 2-colour dark-G instrument artifact (G is the base call for no signal). A C>A or G>T excess
# instead points at oxidative damage during library preparation, and `--min-reads` will not fix it
# -- damage predates the barcode, so every read of the molecule agrees. Enzymatic repair before
# ligation is the fix, and it belongs in the wet lab.
#
# The CLI spelling is `--pre-amp-error`; `--rt-error` is kept as an alias because it shipped.


@dataclass(frozen=True)
class Assay:
    layouts: tuple[str, ...]
    sensitivity: str            # counting | sensitive | ultrasensitive
    min_reads: int
    pre_amp_error: str
    payload_diverse: bool
    extra: tuple[str, ...]      # further flags the assay implies, verbatim
    note: str


ASSAYS: dict[str, Assay] = {
    "airr": Assay(
        ("migec", "primerid"), "counting", 1, "rt", True, (),
        "Repertoire sequencing, and viral quasispecies amplicons like Primer ID, which pose the "
        "same problem: many distinct variants of one locus. The UMI is inline in the read that "
        "also carries the template, and the template is a rearrangement, so the payload is a "
        "second identifier. A clonotype seen once is still a clonotype -- keep --min-reads 1. "
        "Raise it to 3 only when the deliverable is a somatic hypermutation call, not a count.",
    ),
    "amplicon": Assay(
        ("umi",), "sensitive", 2, "rt", False, (),
        "Targeted amplicon panel: a few regions amplified by PCR, read for variants. Same purpose "
        "and same tier as `exome` -- what differs is depth per molecule, because every read lands "
        "on one of a handful of loci rather than being spread over a capture space. Families are "
        "therefore deep, so --min-reads 3 costs little here and is worth taking; below ~1% VAF use "
        "the `ctdna` profile instead, which is what that regime actually needs.",
    ),
    "exome": Assay(
        ("umi",), "sensitive", 2, "rt", False, (),
        "Hybrid-capture exome or gene panel on germline-to-subclonal frequencies. Families are "
        "small (capture duplication is a few-fold, not the hundreds an amplicon gives), so 3 costs "
        "more than it buys; 2 still means every base was seen twice. The vendor sets the UMI "
        "length -- IDT, Twist and Illumina all differ -- so run `migec suggest` on the FASTQ rather "
        "than trusting a preset here.",
    ),
    "ctdna": Assay(
        ("tso500", "duplex"), "ultrasensitive", 3, "7.37e-5", False, (),
        "Circulating tumour DNA: the input is cell-free DNA from a blood draw, the target is the "
        "tumour-derived fraction of it, and that fraction is what the VAF measures. Molecule-"
        "limited by the draw and artifact-limited below ~1%, so both halves of the limit bind at "
        "once. Measured on certified reference material: reliable to 0.25%, and at --min-reads 1 "
        "the dark-G artifact was ADDITIVE to true positives -- the 0.25% arm read 0.79%. Count "
        "molecules per target, never per panel. Note: there is no reverse transcription anywhere "
        "in this assay; the floor is library-prep damage plus the first PCR cycle.",
    ),
    "mrd": Assay(
        ("migec", "duplex"), "ultrasensitive", 3, "rt", True, (),
        "Minimal residual disease: one known clone, tracked to the lowest frequency the input "
        "allows. The IGH rearrangement is the marker, so the payload identifies the molecule as "
        "well as the barcode does -- which is the one advantage this assay has over ctDNA at the "
        "same depth. The limit is set by input DNA: 1e-6 needs ~3e6 informative molecules, and no "
        "consensus threshold substitutes for them. Note: clinical MRD amplifies GENOMIC DNA, so "
        "there is no reverse transcription and the 1e-4 floor is a deliberately conservative "
        "bracket, not an RT rate -- it caps quality lower, which is the safe direction. Drop to "
        "`medium` only if the protocol's polymerase fidelity is known. An RNA-based 5'-RACE MRD "
        "assay does have an RT step, and 1e-4 is then the literal figure.",
    ),
    "rnaseq": Assay(
        ("smarter-umi", "umi"), "counting", 1, "rt", False, ("--fast",),
        "UMI-tagged bulk RNA-seq where the deliverable is a molecule count per gene. --fast is the "
        "counting consensus: modal exact sequence, per-base max quality over the reads carrying "
        "it. It is not a cheaper consensus for variant work -- it is the right one when you are "
        "deduplicating rather than error-correcting.",
    ),
    "10x-gex": Assay(
        ("10x",), "counting", 1, "rt", False, ("--fast",),
        "10x Chromium 3' gene expression. Shallow by design -- 1-3 reads per (cell, UMI) is the "
        "normal case, so a read threshold would delete most of the library. Never: alevin, "
        "bustools and STARsolo must not see a consensus FASTQ; they deduplicate from a raw barcode "
        "read that no longer exists. Feed the consensus to an aligner and count from the BAM.",
    ),
    "10x-vdj": Assay(
        ("10x-v2", "10x"), "counting", 1, "rt", True, ("--contig",),
        "10x Chromium 5' V(D)J. Reads under one (cell, UMI) are random-primed fragments of one "
        "receptor and are NOT co-terminal, so --contig places them and emits one consensus per "
        "overlap component. Never bridged across a gap. Assembling the full-length receptor, "
        "calling doublets and dropping contaminating chains are arda's job, not migec's.",
    ),
}

# What people actually type. Never: "amplicon" is ambiguous in the wild and must not be an alias of
# `airr` -- a targeted panel of a few PCR-amplified regions is also an amplicon assay, and it wants
# the opposite settings (variant calling on a uniform payload, not counting a diverse one). They
# are two profiles, and `amplicon` is the targeted one because that is what the word means outside
# immunology.
ALIASES = {"repseq": "airr", "rep-seq": "airr", "quasispecies": "airr",
           "targeted": "amplicon", "panel": "amplicon", "capture": "exome",
           "10x": "10x-gex", "gex": "10x-gex", "vdj": "10x-vdj",
           "cfdna": "ctdna"}


def assay(name: str) -> tuple[str, Assay]:
    """Look up an assay profile by name or alias. Raises ValueError listing every one."""
    key = ALIASES.get(name.strip().lower(), name.strip().lower())
    if key not in ASSAYS:
        listed = "\n".join(f"  {n:<10} {ASSAYS[n].sensitivity}" for n in ASSAYS)
        raise ValueError(f"unknown assay {name!r}. Available:\n{listed}")
    return key, ASSAYS[key]


def format_assay(name: str) -> str:
    """The paste-ready recipe for one assay, for `migec sheet --assay`."""
    key, a = assay(name)
    layout = a.layouts[0]
    master, slave, _ = PRESETS[layout]
    extra = "".join(" " + f for f in a.extra)
    payload = ("diverse -- a second identifier" if a.payload_diverse
               else "uniform -- the barcode carries the whole burden")
    out = [
        f"{key}  ({a.sensitivity})",
        f"    {a.note}",
        "",
        f"    layout          {layout}   {master}" + (f"   slave {slave}" if slave else ""),
        f"    also fits       {', '.join(a.layouts[1:]) or '-'}",
        f"    payload         {payload}",
        "",
        f"    migec checkout READS.fq.gz --bc-pattern '{master}' --sample S1 -o co/",
        "    migec refine co/S1.fq.gz -o rf/",
        f"    migec assemble rf/S1.fq.gz -o as/ --min-reads {a.min_reads} "
        f"--pre-amp-error {a.pre_amp_error}{extra}",
    ]
    return "\n".join(out)


def format_assays() -> str:
    """Every assay profile, one block each."""
    return "\n\n".join(format_assay(n) for n in ASSAYS)


def format_presets() -> str:
    """The preset table, for `migec sheet --presets` and the docs."""
    used: dict[str, list[str]] = {}
    for name, a in ASSAYS.items():
        for layout in a.layouts:
            used.setdefault(layout, []).append(name)
    out = []
    for name, (master, slave, description) in PRESETS.items():
        out.append(f"{name}")
        out.append(f"    pattern  {master}" + (f"    slave  {slave}" if slave else ""))
        out.append(f"    {description}")
        out.append(f"    assays   {', '.join(used.get(name, [])) or '-'}"
                   f"   (`migec sheet --assay NAME` for the settings each implies)")
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
