"""The migec command line.

Five pipeline commands, and that is the whole surface:

    checkout    extract sample/cell/UMI barcodes, trim, write per-sample FASTQ + QC tables
    suggest     infer where the UMI/primer/cell barcode actually is in the reads
    refine      estimate error rates, correct barcodes, write QC tables and plots
    assemble    build consensus sequences per molecule, write FASTQ
    subsample   take N whole UMIs (with all their reads) to make an example fixture

plus `info` and `sheet`, which read no data and produce no pipeline output.

Adding a sixth pipeline command requires a benchmark the existing five cannot pass. See CLAUDE.md.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from migec import __version__

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="UMI barcode extraction, correction and consensus assembly.",
)


@app.command()
def info() -> None:
    """Print versions and the on-disk format this build reads."""
    from migec import _core

    typer.echo(f"migec           {__version__}")
    typer.echo(f"migec._core     {_core.__version__}")
    typer.echo(f"mig format      v{_core.MIG_FORMAT_VERSION}")


@app.command()
def checkout(
    reads: Path = typer.Argument(..., help="Input FASTQ (R1), optionally gzipped."),
    reads2: Optional[Path] = typer.Argument(
        None,
        help="Second mate. When given, the tag is looked for in either mate and the pair is "
        "swapped so R1 always carries it.",
    ),
    barcodes: Optional[Path] = typer.Option(
        None, "--barcodes", "-b", help="MIGEC-style barcode table. Use this OR --bc-pattern."
    ),
    bc_pattern: Optional[str] = typer.Option(
        None,
        "--bc-pattern",
        help="Where the barcode is, for one sample, inline. Two spellings: a PATTERN, where `N` is "
        "a UMI base and `X` a cell barcode base -- `^NNNNNNNN`, `^NNNNXNNN` -- or a half-open "
        "SLICE list, `0:8` or `0:4,5:10` or `cell:0:16,16:26`. A leading `^`, and any slice list, "
        "anchors the barcode at the first base and needs no --max-offset. Note: umi_tools writes "
        "the cell barcode as `C`, which is cytosine here -- see the error you get if you paste one.",
    ),
    preset: Optional[str] = typer.Option(
        None,
        "--preset",
        help="A named layout instead of a pattern: umi, migec, primerid, duplex, 10x, 10x-v2, "
        "tso500, smarter-umi. `migec sheet --presets` prints each one with what it is and where "
        "the layout is written down.",
    ),
    read_structure: Optional[str] = typer.Option(
        None,
        "--read-structure",
        help="An fgbio/Picard read structure instead of a pattern, which is what TSO500, fgbio "
        "and samtools all speak. M is a UMI base, B a sample/cell barcode, S skipped, T template: "
        "TSO500 is 5M5S+T and 10x 5' is 16B10M+T. Positional by definition, so it carries its own anchor.",
    ),
    read_structure2: Optional[str] = typer.Option(
        None,
        "--read-structure2",
        help="Read structure for the second mate, when it carries barcode too. TSO500 puts a UMI "
        "on both mates, and the two halves concatenate into one molecule identifier.",
    ),
    sample_id: str = typer.Option(
        "sample", "--sample", help="Sample name for --bc-pattern. Ignored with --barcodes."
    ),
    out_dir: Path = typer.Option(..., "--out", "-o", help="Output directory."),
    threads: int = typer.Option(
        0,
        "--threads",
        "-t",
        help="Worker threads; 0 uses one per core. Output is byte-identical whatever this is.",
    ),
    trim: str = typer.Option(
        "pattern",
        "--trim",
        help="'pattern' drops the adapter, sample tag and UMI, leaving the payload; "
        "'none' keeps the read whole and puts the UMI in the header only.",
    ),
    min_umi_quality: int = typer.Option(
        0,
        "--min-umi-quality",
        help="Drop reads whose worst UMI base is below this Phred. 0 (default) keeps everything: "
        "a low-quality UMI base is a reason to be less certain, not to discard the molecule.",
    ),
    max_offset: Optional[int] = typer.Option(
        None,
        "--max-offset",
        help="Where the pattern may start. Default is automatic: a layout with nothing to score "
        "is anchored at the first base (0), anything with an adapter to place it gets a free scan "
        "(-1). Set it only to override that. A short anchor like a 5 nt dual-end handle occurs by "
        "chance every kilobase, so a free scan cannot place it and correctly refuses to.",
    ),
    write_unmatched: bool = typer.Option(
        False, "--write-unmatched", help="Also write reads that matched no pattern."
    ),
    limit_read: int = typer.Option(
        0,
        "--limit-read",
        help="Stop after this many input reads. A smoke test on a big run -- never a sample: the "
        "first N reads of a FASTQ are one corner of one flowcell, so nothing measured under a "
        "limit describes the library. Use `migec subsample` when you want a fixture.",
    ),
    mig: bool = typer.Option(
        False,
        "--mig",
        help="Write `<sample>.<bucket>.mig` buckets instead of one FASTQ per sample. The reads "
        "come out already range-partitioned on the barcode, which is the pass `migec assemble` "
        "otherwise spends most of its time in -- point assemble at the output directory and it "
        "skips it. A `.mig` file is a migec intermediate: nothing else reads it, so keep the "
        "default when the reads are going to an aligner.",
    ),
) -> None:
    """Demultiplex by barcode pattern, extract and trim UMIs, write QC tables."""
    from migec.checkout import format_report, run

    if trim not in ("pattern", "none"):
        raise typer.BadParameter("--trim must be 'pattern' or 'none'")
    from migec.sheet import from_read_structure, parse_layout

    # Checked first: the second mate's structure alone is a truncated command, and saying "give
    # exactly one of four things" to someone who gave one of them is an unhelpful place to land.
    if read_structure2 is not None and read_structure is None:
        raise typer.BadParameter("--read-structure2 needs --read-structure")
    given = [n for n, v in (("--barcodes", barcodes), ("--bc-pattern", bc_pattern),
                            ("--read-structure", read_structure), ("--preset", preset)) if v]
    if len(given) != 1:
        raise typer.BadParameter(
            "give exactly one of --barcodes, --bc-pattern, --read-structure and --preset"
            + (f" (got {', '.join(given)})" if given else "")
        )

    slave, anchored = None, False
    try:
        if preset is not None:
            from migec.sheet import preset as lookup

            bc_pattern, slave_spec = lookup(preset)
            if slave_spec:
                slave, anchored = parse_layout(slave_spec)
        if read_structure is not None:
            # A read structure is positional by definition -- every count is measured from the
            # first base -- so it carries its own anchor and never needs --max-offset spelled out.
            bc_pattern = from_read_structure(read_structure)
            slave = from_read_structure(read_structure2) if read_structure2 else None
            anchored = True
        if bc_pattern is not None:
            bc_pattern, master_anchored = parse_layout(bc_pattern)
            anchored = anchored or master_anchored
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc

    # A caret, a slice list or a read structure all say the same thing: the barcode is AT a
    # position, so there is nothing to scan for. Spelling `--max-offset 0` after that is
    # boilerplate the user cannot get right without knowing why it exists.
    if max_offset is None and anchored:
        max_offset = 0
    if bc_pattern is not None:
        barcodes = _inline_sheet(bc_pattern, sample_id, out_dir, slave)
    summary = run(
        reads,
        barcodes,
        out_dir,
        reads2=reads2,
        trim=trim,
        min_umi_quality=min_umi_quality,
        write_unmatched=write_unmatched,
        threads=threads,
        max_offset=max_offset,
        limit_reads=limit_read,
        mig=mig,
    )
    typer.echo(format_report(summary))
    if mig:
        typer.echo(
            f"\nwrote {len(summary['mig_paths'])} .mig bucket(s) in {out_dir}"
            f"\n      assemble a sample with `migec assemble {out_dir}/"
            f"{summary['samples'][0]['sample_id']}.000.mig -o <out>` -- one bucket names the "
            f"whole partition, and the partition pass is already done"
        )
    typer.echo(
        f"\nwrote {out_dir}/checkout.{{summary,coverage,umi_composition,barcode_space,"
        f"umi_quality,quality_calibration,pattern_positions,trimming}}.tsv"
        f"\n      draw them with `migec plot {out_dir}`"
    )


@app.command()
def sheet(
    barcodes: Optional[Path] = typer.Argument(None, help="MIGEC-style barcode table to inspect."),
    presets: bool = typer.Option(
        False, "--presets", help="List the named layouts instead, and where each one is from."
    ),
    assay: Optional[str] = typer.Option(
        None,
        "--assay",
        help="Print the paste-ready recipe for an assay: amplicon (airr), exome, ctdna, mrd, "
        "rnaseq, 10x-gex, 10x-vdj. Pass `all` for every one. A layout says where the barcode is; "
        "an assay also says what a consensus is worth.",
    ),
) -> None:
    """Show what each row of a barcode table will extract, without running anything."""
    from migec.sheet import describe, format_assay, format_assays, format_presets, read_barcodes

    if assay:
        try:
            typer.echo(format_assays() if assay.strip().lower() == "all" else format_assay(assay))
        except ValueError as exc:
            raise typer.BadParameter(str(exc)) from exc
        return
    if presets:
        typer.echo(format_presets())
        return
    if barcodes is None:
        raise typer.BadParameter(
            "give a barcode table, --presets to list the named layouts, or --assay NAME "
            "for what an experiment implies"
        )
    typer.echo(describe(read_barcodes(barcodes)))


@app.command()
def plot(
    directory: Path = typer.Argument(
        ..., help="A stage's output directory. Every table it finds there gets a figure."
    ),
    out_dir: Optional[Path] = typer.Option(
        None, "--out", "-o", help="Where to write. Defaults to <directory>/plots."
    ),
    fmt: str = typer.Option("svg", "--format", help="svg, png or pdf."),
) -> None:
    """Draw the QC figures from the tables a stage already wrote.

    Reads no reads and produces no pipeline output: every panel is a gnuplot script over a
    committed TSV, so a figure can be redrawn from the table next to it long after the FASTQ is
    gone. gnuplot itself is not a Python package -- without it the scripts are still written.
    """
    from migec.plot import format_report, run

    try:
        summary = run(directory, out_dir, fmt=fmt)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(format_report(summary))


@app.command()
def suggest(
    reads: Path = typer.Argument(..., help="FASTQ to profile. For paired data try both mates."),
    out_dir: Optional[Path] = typer.Option(
        None, "--out", "-o", help="Write suggest.cycles.tsv, suggest.segments.tsv, suggest.json."
    ),
    cycles: int = typer.Option(60, "--cycles", help="Leading cycles to profile."),
    max_reads: int = typer.Option(
        200_000, "--max-reads", help="Reads to profile. The composition converges long before this."
    ),
    umi_deviation: float = typer.Option(
        0.18,
        "--umi-deviation",
        help="How far a cycle may sit from a flat 1/4 base composition and still be called UMI. "
        "Real synthesiser mixes are routinely 20/30/30/20, which is 0.05.",
    ),
) -> None:
    """Infer where the UMI and primer are, from the per-cycle base composition."""
    from migec.suggest import format_report, run

    summary = run(reads, out_dir, cycles=cycles, max_reads=max_reads, umi_deviation=umi_deviation)
    typer.echo(format_report(summary))


@app.command()
def refine(
    reads: Path = typer.Argument(
        ...,
        help="A per-sample FASTQ written by `migec checkout` -- or the `.mig` buckets that "
        "`checkout --mig` wrote, given as a directory or as one bucket file, in which case refine "
        "writes buckets back and `migec assemble` takes them straight from here.",
    ),
    out_dir: Path = typer.Option(..., "--out", "-o", help="Output directory."),
    sample_id: str = typer.Option("", "--sample", help="Defaults to the BC tag in the reads."),
    min_posterior: float = typer.Option(
        0.95,
        "--min-posterior",
        help="Posterior above which a barcode is folded into its parent. Raise it to correct "
        "less: a wrong merge deletes a molecule and nothing downstream can tell, while a missed "
        "correction only inflates a count.",
    ),
    expect_cells: int = typer.Option(
        3000,
        "--expect-cells",
        help="Cells expected, for the OrdMag call: the 99th percentile of the top this many "
        "barcodes, over ten. Ignored without cell barcodes. EmptyDrops-style rescue of low-count "
        "cells is Cell Ranger's job and is deliberately not reproduced.",
    ),
    cell_whitelist: Optional[Path] = typer.Option(
        None,
        "--cell-whitelist",
        help="List of the cell barcodes that were actually synthesised (10x ships one). Observed "
        "barcodes one substitution away are snapped to it, weighed against the measured prior "
        "that a barcode is genuinely off-list -- without which every hopped or undeclared barcode "
        "is absorbed into whichever entry happens to be nearest.",
    ),
    target_fdr: float = typer.Option(
        0.05,
        "--target-fdr",
        help="Residual false-molecule rate the reported MIG size threshold aims at. The threshold "
        "is REPORTED, never applied: a molecule seen three times with no plausible parent is "
        "information, and cutting it discards real sequence.",
    ),
    threads: int = typer.Option(
        0, "--threads", "-t",
        help="Worker threads for the neighbourhood scan, the error-rate estimate, the residual-FDR "
        "scan and the read rewrite; 0 uses one per core. The output is byte-identical whatever "
        "this is: all four are pure functions of the barcode table, and the merges the scan finds "
        "are applied serially afterwards.",
    ),
    limit_read: int = typer.Option(
        0,
        "--limit-read",
        help="Stop after this many input reads. A smoke test on a big run -- never a sample: the "
        "first N reads of a FASTQ are one corner of one flowcell, so nothing measured under a "
        "limit describes the library. Use `migec subsample` when you want a fixture.",
    ),
    limit_umi: int = typer.Option(
        0,
        "--limit-umi",
        help="Stop once this many distinct barcodes have been seen. Same warning as --limit-read: "
        "these are the barcodes that happen to appear first, not a sample of them.",
    ),
    no_quality: bool = typer.Option(
        False, "--no-quality", help="Ignore the barcode's own base quality (QX)."
    ),
    no_payload: bool = typer.Option(
        False,
        "--no-payload",
        help="Ignore payload agreement. Both flags exist to measure what the count ratio alone "
        "would have done, which is all the first version had.",
    ),
) -> None:
    """Correct barcode errors and rewrite the reads with the corrected barcode."""
    from migec.refine import format_report, run

    summary = run(
        reads, out_dir, sample_id=sample_id, use_quality=not no_quality,
        use_payload=not no_payload, min_posterior=min_posterior, expect_cells=expect_cells,
        cell_whitelist=cell_whitelist or "", target_fdr=target_fdr, threads=threads,
        limit_reads=limit_read, limit_umis=limit_umi,
    )
    typer.echo(format_report(summary))
    if summary.get("mig_paths"):
        typer.echo(
            f"\nwrote {len(summary['mig_paths'])} .mig bucket(s), re-partitioned on the "
            f"corrected barcode"
            f"\n      assemble them with `migec assemble {out_dir}/{summary['sample_id']}.000.mig "
            f"-o <out>`"
        )


@app.command()
def assemble(
    reads: Path = typer.Argument(
        ...,
        help="A per-sample FASTQ written by `migec checkout` -- or the `.mig` buckets that "
        "`checkout --mig` wrote, given as a directory or as one bucket file, in which case the "
        "partition pass is skipped because it has already been done.",
    ),
    out_dir: Path = typer.Option(..., "--out", "-o", help="Output directory."),
    sample_id: str = typer.Option("", "--sample", help="Defaults to the BC tag in the reads."),
    rt_error: str = typer.Option(
        "rt",
        "--pre-amp-error",
        "--rt-error",
        help="The pre-amplification error floor, which caps every emitted quality at -10 log10 of "
        "it. Anything already in the molecule when the first amplification cycle started is in "
        "every read of the group, so no consensus removes it. Name the chemistry: 'rt' (1e-4, "
        "Q40), 'medium' (1e-5, Q50) for an ordinary polymerase, 'high' (1e-6, Q60) for a "
        "proofreading one -- or give the rate itself, or 'auto' to FIT it from this dataset's own "
        "deep molecules. 'auto' needs a clonal template and costs a second assembly pass, and it "
        "says which class it fell back to when the library cannot support the fit. Note: the "
        "class names are historical. Only "
        "an RNA assay has a reverse transcription step; on a DNA library the same floor is set by "
        "library-prep damage (oxidation during shearing, cytosine deamination) plus the first PCR "
        "cycle, and 'rt' is then just the name of the 1e-4 bracket. It is the ONE-MOLECULE floor "
        "and every record here is one molecule: 10x's Q60 needs two UMIs to agree, and combining "
        "molecules is arda's job. See docs/quality_floor.rst for what measured each.",
    ),
    contig: bool = typer.Option(
        False,
        "--contig",
        help="Random-primed reads sharing a barcode tile the molecule instead of starting at the "
        "same base. Place them against each other and emit one consensus per overlap component. "
        "X1 measured 27.3% of 10x groups holding more than one component (docs/fragmented.rst); "
        "one consensus over those asserts sequence no read covers.",
    ),
    fast: bool = typer.Option(
        False,
        "--fast",
        help="Counting mode: emit each group's most frequent EXACT sequence, with every base "
        "carrying the best quality any read of that sequence reported. No column model, so no "
        "per-base error correction and no sub-clustering -- use it when the deliverable is "
        "molecule COUNTS (expression, clonotype abundance) rather than error-free sequence. "
        "Incompatible with --contig, whose reads tile the molecule and share no exact sequence.",
    ),
    threads: int = typer.Option(
        0, "--threads", "-t",
        help="Workers for the consensus pass; 0 uses one per core. Buckets are the unit of work "
        "and there are always at least 16 of them, so the output is byte-identical whatever this "
        "is -- bucket order is barcode order, and the bucket count does not depend on -t.",
    ),
    limit_read: int = typer.Option(
        0,
        "--limit-read",
        help="Stop after this many input reads. A smoke test on a big run -- never a sample: the "
        "first N reads of a FASTQ are one corner of one flowcell, so nothing measured under a "
        "limit describes the library. Use `migec subsample` when you want a fixture.",
    ),
    limit_umi: int = typer.Option(
        0,
        "--limit-umi",
        help="Stop once this many distinct barcodes have been seen. Same warning as --limit-read: "
        "these are the barcodes that happen to appear first, not a sample of them.",
    ),
    min_reads: int = typer.Option(
        1,
        "--min-reads",
        help="Molecules below this are dropped. The default keeps everything: a molecule seen "
        "three times is still sequence.",
    ),
) -> None:
    """Collapse the reads of each barcode (sample + cell + UMI) into a consensus."""
    from migec.assemble import format_report, parse_rt_error, run

    try:
        rt_floor = parse_rt_error(rt_error)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    # Refused here, on the caller's thread, rather than producing one consensus per read: two
    # fragments of one molecule are different strings by construction, so a modal vote over them
    # returns whichever fragment was seen most and silently drops the rest of the molecule.
    if fast and contig:
        raise typer.BadParameter(
            "--fast and --contig are incompatible. --contig places reads that TILE a molecule, "
            "and tiling reads share no exact sequence for --fast to take a majority over"
        )
    summary = run(
        reads, out_dir, sample_id=sample_id, rt_floor=rt_floor, contig=contig, fast=fast,
        min_reads=min_reads, threads=threads, limit_reads=limit_read, limit_umis=limit_umi,
    )
    typer.echo(format_report(summary))


@app.command()
def subsample(
    reads: Path = typer.Argument(..., help="A FASTQ carrying RX (and CB) tags."),
    output: Path = typer.Option(..., "--out", "-o", help="Output FASTQ; .gz is honoured."),
    keep: float = typer.Option(
        1.0,
        "--keep",
        help="Percent of BARCODES to keep, with all of their reads. Never a percent of the reads: "
        "sampling reads gives one read per molecule and destroys the MIG size distribution, which "
        "is the one thing a UMI fixture exists to preserve.",
    ),
    by_umi_only: bool = typer.Option(
        False,
        "--by-umi-only",
        help="Sample molecules rather than whole cells. Off by default, because a fixture of "
        "thousands of cells holding one molecule each is the same mistake as sampling reads.",
    ),
) -> None:
    """Take all the reads of a fraction of the barcodes."""
    from migec.subsample import format_report, run

    summary = run(reads, output, keep_percent=keep, by_cell=not by_umi_only)
    typer.echo(format_report(summary))


def _inline_sheet(
    pattern: str, sample_id: str, out_dir: Path, slave: Optional[str] = None
) -> Path:
    """Write a one-row barcode table for `--bc-pattern`, after refusing an ambiguous one.

    Never: umi_tools spells a cell-barcode position `C`, and `C` is cytosine in this grammar. Pasting
    a umi_tools pattern would compile -- into a pattern demanding a run of literal cytosines,
    which matches nothing and looks like a bad library rather than a bad flag. Refuse it and say
    what to write instead.
    """
    import re

    run_of_c = max((len(m) for m in re.findall(r"C+", pattern)), default=0)
    if run_of_c >= 4 and "X" not in pattern.upper():
        raise typer.BadParameter(
            f"--bc-pattern {pattern!r} contains a run of {run_of_c} literal C. In this grammar C "
            f"is cytosine, not a cell barcode -- umi_tools spells it that way, migec spells it X. "
            f"Write {pattern.replace('C', 'X')!r} if you meant a cell barcode, or keep the C's "
            f"and add --max-offset if you really do mean a run of cytosines."
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "bc_pattern.txt"
    row = f"{sample_id}\t{pattern}" + (f"\t{slave}" if slave else "")
    path.write_text(row + "\n")
    return path
