"""The migec command line.

Five commands, and that is the whole surface:

    checkout    extract sample/cell/UMI barcodes, write .mig buckets
    suggest     infer where the UMI/primer/cell barcode actually is in the reads
    refine      estimate error rates, correct barcodes, write QC tables and plots
    assemble    build consensus sequences per molecule, write FASTQ
    subsample   take N whole UMIs (with all their reads) to make an example fixture

Adding a sixth requires a benchmark the existing five cannot pass. See CLAUDE.md.
"""

from __future__ import annotations

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
def checkout() -> None:
    """Extract barcodes from FASTQ into .mig buckets. (M2)"""
    raise typer.Exit(_not_yet("checkout", "M2"))


@app.command()
def suggest() -> None:
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
        f"`migec {name}` is not implemented yet (planned for {milestone}).\n"
        f"This build ships the .mig format and FASTQ IO; see ROADMAP.md.",
        err=True,
    )
    return 2
