"""Regenerate the example QC tables and figures the README embeds.

    python scripts/example_figures.py            # writes assets/*.tsv and assets/*.svg

2026-08-13. The four panels the README shows -- the barcode rank plot, the MIG size spectrum, the
rank/Zipf curve and consensus quality against depth -- need a library with the right *shape*, not a
big one: a real cell/ambient split so the knee exists, a heavy-tailed amplification so the Zipf
curve has slope, and per-base errors so the emitted quality has spread. A uniform simulator draws
all four as straight lines and the figures say nothing.

The tables it writes are committed; the SVGs are drawn from those tables by `migec plot`, so a
figure can always be redrawn from the numbers next to it. Provenance is in SOURCES.md: these are
DERIVED (simulated), never experimental.
"""

from __future__ import annotations

import argparse
import gzip
import random
import shutil
import tempfile
from pathlib import Path

from migec.assemble import run as assemble_run
from migec.checkout import run as checkout_run
from migec.plot import run as plot_run
from migec.refine import run as refine_run

SEED = 3
CELLS, AMBIENT = 120, 4000
# A spread of reported Phred rather than a constant, so the consensus quality panel has quartiles
# to draw. Roughly what a 2-colour instrument emits: a few distinct values, most of them high.
QUALS = "5:?ABFHI"
# Pareto alpha 1.1: heavy enough that the rank curve has a visible slope over three decades. A
# uniform or Poisson depth gives a flat line and the panel is then a test of nothing.
PARETO_ALPHA = 1.1
CELL_LEN, UMI_LEN, PAYLOAD_LEN = 16, 10, 60


def simulate(directory: Path, rng: random.Random) -> tuple[Path, Path, Path]:
    def mutate(seq: str, rate: float) -> str:
        return "".join(rng.choice("ACGT") if rng.random() < rate else b for b in seq)

    def barcode(n: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(n))

    r1, r2, sheet = directory / "r1.fq.gz", directory / "r2.fq.gz", directory / "bc.txt"
    n = 0
    with gzip.open(r1, "wt") as f1, gzip.open(r2, "wt") as f2:
        def emit(cell: str, umi: str, payload: str, reads: int) -> None:
            nonlocal n
            for _ in range(reads):
                quality = "".join(rng.choice(QUALS) for _ in range(PAYLOAD_LEN))
                f1.write(f"@r{n}\n{mutate(cell + umi, 0.002)}\n+\n{'I' * 26}\n")
                f2.write(f"@r{n}\n{mutate(payload, 0.01)}\n+\n{quality}\n")
                n += 1

        for _ in range(CELLS):  # real cells: tens to hundreds of molecules, amplified unevenly
            cell = barcode(CELL_LEN)
            for _ in range(rng.randint(40, 200)):
                emit(cell, barcode(UMI_LEN), barcode(PAYLOAD_LEN),
                     max(1, int(rng.paretovariate(PARETO_ALPHA))))
        for _ in range(AMBIENT):  # ambient: a molecule or two, seen once. This makes the knee.
            cell = barcode(CELL_LEN)
            for _ in range(rng.randint(1, 3)):
                emit(cell, barcode(UMI_LEN), barcode(PAYLOAD_LEN), 1)

    sheet.write_text("PBMC\t" + "X" * CELL_LEN + "N" * UMI_LEN + "\n")
    return r1, r2, sheet


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--out", default="assets", type=Path)
    parser.add_argument("--keep", action="store_true", help="keep the simulated FASTQ")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    work = Path(tempfile.mkdtemp(prefix="migec-example-"))
    try:
        rng = random.Random(SEED)
        r1, r2, sheet = simulate(work, rng)
        checkout_run(r1, sheet, work / "co", reads2=r2)
        refine_run(work / "co" / "PBMC_R2.fq.gz", work / "ref", sample_id="PBMC",
                   expect_cells=CELLS)
        assemble_run(work / "ref" / "PBMC.fq.gz", work / "asm", sample_id="PBMC")

        # Only the tables the README's panels need. Never the whole output directory: a table with
        # no figure is a file nobody will ever notice going stale.
        # Never PBMC.mig.tsv: assets/ already holds SRR1763769.mig.tsv, which is REAL data, and
        # two matches for the same glob draw the per-molecule panels twice -- once from a
        # simulation and once from a library, side by side with nothing saying which is which.
        for name in ("ref/PBMC.cell_rank.tsv", "ref/PBMC.sizes.tsv", "ref/PBMC.umi_errors.tsv",
                     "asm/assemble.quality_by_depth.tsv", "asm/assemble.coverage.tsv"):
            shutil.copy(work / name, args.out / Path(name).name)
        if args.keep:
            shutil.copy(r1, args.out / "example_r1.fq.gz")
            shutil.copy(r2, args.out / "example_r2.fq.gz")
    finally:
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)

    summary = plot_run(args.out, args.out)
    print(f"tables and {len(summary['drawn'])} figures in {args.out}")
    for name in summary["drawn"]:
        print(f"  {name}")
    for failure in summary["failed"]:
        print(f"warning: {failure}")


if __name__ == "__main__":
    main()
