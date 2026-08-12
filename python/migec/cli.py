"""The migec command line.

Five commands, and that is the whole surface:

    checkout    extract sample/cell/UMI barcodes, trim, write per-sample FASTQ + QC tables
    suggest     infer where the UMI/primer/cell barcode actually is in the reads
    refine      estimate error rates, correct barcodes, write QC tables and plots
    assemble    build consensus sequences per molecule, write FASTQ
    subsample   take N whole UMIs (with all their reads) to make an example fixture

Adding a sixth requires a benchmark the existing five cannot pass. See CLAUDE.md.
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
    reads: Optional[Path] = typer.Argument(None),
) -> None:
    """Infer UMI/primer/cell-barcode placement from the reads. (M4)"""
    raise typer.Exit(_not_yet("suggest", "M4"))


@app.command()
def refine() -> None:
    """Estimate error rates, correct barcodes, emit QC tables. (M3)"""
    raise typer.Exit(_not_yet("refine", "M3"))


@app.command()
def assemble() -> None:
    """Assemble consensus sequences per molecule. (M1)"""
    raise typer.Exit(_not_yet("assemble", "M1"))


@app.command()
def subsample() -> None:
    """Take N whole UMIs with all of their reads. (M4)"""
    raise typer.Exit(_not_yet("subsample", "M4"))


def _not_yet(name: str, milestone: str) -> int:
    typer.echo(
        f"`migec {name}` is not implemented yet (planned for {milestone}). See ROADMAP.md.",
        err=True,
    )
    return 2
