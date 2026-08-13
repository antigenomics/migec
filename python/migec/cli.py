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
    )
    typer.echo(format_report(summary))
    typer.echo(f"\nwrote {out_dir}/checkout.{{summary,coverage,umi_composition}}.tsv")


@app.command()
def sheet(
    barcodes: Optional[Path] = typer.Argument(None, help="MIGEC-style barcode table to inspect."),
    presets: bool = typer.Option(
        False, "--presets", help="List the named layouts instead, and where each one is from."
    ),
) -> None:
    """Show what each row of a barcode table will extract, without running anything."""
    from migec.sheet import describe, format_presets, read_barcodes

    if presets:
        typer.echo(format_presets())
        return
    if barcodes is None:
        raise typer.BadParameter("give a barcode table, or --presets to list the named layouts")
    typer.echo(describe(read_barcodes(barcodes)))


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
    reads: Path = typer.Argument(..., help="A per-sample FASTQ written by `migec checkout`."),
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
        cell_whitelist=cell_whitelist or "", target_fdr=target_fdr,
    )
    typer.echo(format_report(summary))


@app.command()
def assemble(
    reads: Path = typer.Argument(..., help="A per-sample FASTQ written by `migec checkout`."),
    out_dir: Path = typer.Option(..., "--out", "-o", help="Output directory."),
    sample_id: str = typer.Option("", "--sample", help="Defaults to the BC tag in the reads."),
    rt_error: float = typer.Option(
        1e-4,
        "--rt-error",
        help="The RT/first-cycle-PCR error floor, which caps every emitted quality at "
        "-10 log10 of it. Default measured on an HIV-1 Primer ID control (docs/quality_floor.rst); "
        "the 1e-6 that gets assumed is excluded by two orders of magnitude.",
    ),
    contig: bool = typer.Option(
        False,
        "--contig",
        help="Random-primed reads sharing a barcode tile the molecule instead of starting at the "
        "same base. Place them against each other and emit one consensus per overlap component. "
        "X1 measured 27.3% of 10x groups holding more than one component (docs/fragmented.rst); "
        "one consensus over those asserts sequence no read covers.",
    ),
    min_reads: int = typer.Option(
        1,
        "--min-reads",
        help="Molecules below this are dropped. The default keeps everything: a molecule seen "
        "three times is still sequence.",
    ),
) -> None:
    """Collapse the reads of each barcode (sample + cell + UMI) into a consensus."""
    from migec.assemble import format_report, run

    summary = run(
        reads, out_dir, sample_id=sample_id, rt_floor=rt_error, contig=contig,
        min_reads=min_reads,
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


def _not_yet(name: str, milestone: str) -> int:
    typer.echo(
        f"`migec {name}` is not implemented yet (planned for {milestone}). See ROADMAP.md.",
        err=True,
    )
    return 2
