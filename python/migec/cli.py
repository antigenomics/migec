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
    barcodes: Path = typer.Option(..., "--barcodes", "-b", help="MIGEC-style barcode table."),
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
    max_offset: int = typer.Option(
        -1,
        "--max-offset",
        help="Where the pattern may start. -1 scans the whole read; 0 anchors it at the first "
        "base. Positional chemistries need 0 -- a short anchor like a 5 nt dual-end handle occurs "
        "by chance every kilobase, so a free scan cannot place it and correctly refuses to.",
    ),
    write_unmatched: bool = typer.Option(
        False, "--write-unmatched", help="Also write reads that matched no pattern."
    ),
) -> None:
    """Demultiplex by barcode pattern, extract and trim UMIs, write QC tables."""
    from migec.checkout import format_report, run

    if trim not in ("pattern", "none"):
        raise typer.BadParameter("--trim must be 'pattern' or 'none'")
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
    barcodes: Path = typer.Argument(..., help="MIGEC-style barcode table to inspect."),
) -> None:
    """Show what each row of a barcode table will extract, without running anything."""
    from migec.sheet import describe, read_barcodes

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


def _not_yet(name: str, milestone: str) -> int:
    typer.echo(
        f"`migec {name}` is not implemented yet (planned for {milestone}). See ROADMAP.md.",
        err=True,
    )
    return 2
