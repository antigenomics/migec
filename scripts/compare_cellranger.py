#!/usr/bin/env python3
# 2026-08-14
#
# migec against Cell Ranger 5.0.0 on a 10x droplet library, `sc5p_v2_hs_PBMC_1k` VDJ-T.
#
# The overlap between the two is the droplet front end: given R1 = 16 nt cell barcode + 10 nt UMI,
# which reads carry a real GEM barcode, which barcodes are cells, and what share of the library is
# in them. Three axes, all of which both tools answer for the whole library rather than for a
# receptor. Cell Ranger's own published output for this exact library is the comparator.
#
# Note: the two tools do NOT answer the same question about a cell. migec's gate is molecules of
# any sequence; Cell Ranger's is "assembled a productive V(D)J contig", so a B cell or a monocyte
# is correctly a migec cell and correctly not a Cell Ranger one. The cell sets are therefore
# reported as five counts and never as a ratio -- a bare 479-against-N reads as an accuracy claim
# and is a population statement. The number that says whether the difference costs anything is
# reads in called cells, which is on the same denominator for both.
#
# Never: `Cell Ranger's AIRR export means the opposite of arda's by the same column name.` In
# `_airr_rearrangement.tsv`, `consensus_count` is READS and `duplicate_count` is UMIs -- verified
# on all 943 rows against the contig CSV, and the inverse of what arda writes
# (`arda/src/arda/rnaseq/correct.py:1332`). This script therefore reads counts only from
# `_filtered_contig_annotations.csv`, whose `reads` and `umis` are unambiguous, and asserts the
# identity at runtime rather than trusting it.
#
# Never: `the molecule counts are not comparable and are not published as if they were.` Cell
# Ranger's `umis` counts UMIs incorporated into a filtered contig -- 7,623 over 479 cells -- while
# migec counts every molecule of any sequence, a median of 178 per cell on the same barcodes.
# That is a different population, not an over-count, so no ratio of the two appears here.
#
# Note: `the cost row comes from a run, not from a table.` `metrics_summary.csv` carries no runtime
# or memory, so the Cell Ranger figures are read off `/usr/bin/time -v` in
# `scripts/cellranger_vdj.sbatch` -- Cell Ranger is Linux x86_64, so it runs on the cluster and
# migec's own row is measured on the same input. Pass them with `--cellranger-seconds` and
# `--cellranger-rss-bytes`; omit them and the columns are blank rather than guessed.
#
# Usage:
#     python scripts/compare_cellranger.py \
#         --r1 /tmp/cr/in/R1.fq.gz --r2 /tmp/cr/in/R2.fq.gz \
#         --cellranger-dir ~/data/10x/sc5p_v2_hs_PBMC_1k_t_cellranger \
#         --whitelist ~/data/10x/737K-august-2016.txt \
#         --out /tmp/cr/run --threads 8 --tsv assets/cellranger.tsv
#
# Get the Cell Ranger side (no binary needed -- these are 10x's own published outputs):
#     B=https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k_t
#     for f in filtered_contig_annotations.csv airr_rearrangement.tsv metrics_summary.csv; do
#         curl -sO ${B}_${f}; done
# Get the whitelist (10x ship it inside the cellranger tarball; this is their own mirror of it):
#     curl -L -o 737K-august-2016.txt \
#       https://github.com/10XGenomics/supernova/raw/master/tenkit/lib/python/tenkit/barcodes/737K-august-2016.txt

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Cell Ranger reports this library as 6,301,573 read pairs, which is BOTH lanes. L001 alone is
# 3,155,166 and is what every earlier migec measurement of this library saw.
CELLRANGER_READ_PAIRS = 6_301_573
CELLRANGER_CELLS = 479
CELLRANGER_CONTIGS = 943


def read_whitelist(path: Path) -> set[str]:
    """The 10x barcode list, one per line. Everything from the first `-` is dropped, as
    `migec`'s own loader does (src/whitelist.cpp:22-24) -- 10x write a `-1` gem-group suffix."""
    wl = set()
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        wl.add(line.split("-", 1)[0])
    if not wl:
        raise SystemExit(f"{path} holds no barcodes")
    lengths = {len(b) for b in wl}
    if len(lengths) != 1:
        raise SystemExit(f"{path} mixes barcode lengths {sorted(lengths)} -- it is not a whitelist")
    return wl


def read_cellranger_cells(contig_csv: Path) -> tuple[set[str], int]:
    """The Cell Ranger cell set, and the contig count.

    Never: filter on `is_cell`. `_all_contig_annotations.csv` spans 699 barcodes of which only 479
    are cells, and nothing else in the file says you took the wrong ones.
    """
    cells, contigs = set(), 0
    with open(contig_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["is_cell"].strip().lower() != "true":
                continue
            contigs += 1
            cells.add(row["barcode"].split("-", 1)[0])
    return cells, contigs


def check_airr_counts(airr_tsv: Path, contig_csv: Path) -> None:
    """Read back what Cell Ranger's AIRR count columns actually mean, and refuse if they moved.

    `consensus_count` must equal the contig CSV's `reads` and `duplicate_count` its `umis`. arda
    spells these two the other way round, so a join on the column names silently swaps reads for
    molecules -- a factor of 526 on this library.
    """
    reads, umis = {}, {}
    with open(contig_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            reads[row["contig_id"]] = int(row["reads"])
            umis[row["contig_id"]] = int(row["umis"])
    with open(airr_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = row["sequence_id"]
            if cid not in reads:
                continue
            cc, dc = int(float(row["consensus_count"])), int(float(row["duplicate_count"]))
            if cc != reads[cid] or dc != umis[cid]:
                raise SystemExit(
                    f"Cell Ranger AIRR row {cid}: consensus_count {cc} / duplicate_count {dc} "
                    f"against reads {reads[cid]} / umis {umis[cid]} in the contig CSV. The count "
                    f"columns have moved; arda spells these two the other way round, so nothing "
                    f"downstream of here can be trusted to be counting the same thing"
                )


def read_metrics(path: Path) -> dict[str, str]:
    """`metrics_summary.csv` is two lines, header and values, with quoted thousands separators."""
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    return dict(zip(rows[0], rows[1]))


def _pct(value: str) -> float:
    return float(value.strip().rstrip("%")) / 100.0


def per_cell_reads(barcodes_tsv: Path) -> tuple[dict[str, int], int, int]:
    """Reads per cell, total reads and total molecules, from refine's own barcode table.

    Never: this cannot come from `<sample>.cells.tsv` or `<sample>.cell_rank.tsv` -- every column
    in both is MOLECULES. Ambient barcodes are molecule-rich and read-poor, so the molecule
    fraction reads 57% where the read fraction, which is what Cell Ranger's 86.8% means, reads 84%.
    `<sample>.barcodes.tsv` is one row per (cell, umi) and carries `reads`.
    """
    reads: dict[str, int] = {}
    total_reads = molecules = 0
    with open(barcodes_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n = int(row["reads"])
            reads[row["cell"]] = reads.get(row["cell"], 0) + n
            total_reads += n
            molecules += 1
    return reads, total_reads, molecules


def called_cells(cells_tsv: Path) -> set[str]:
    with open(cells_tsv, newline="") as fh:
        return {r["cell"] for r in csv.DictReader(fh, delimiter="\t") if r["called"] == "1"}


def run_checkout(r1: Path, r2: Path, out: Path, threads: int) -> tuple[dict, Path]:
    """Demultiplex once. The whitelist is a `refine` knob, so checkout does not depend on it and
    is not re-run per whitelist setting -- that was half the wall clock of this comparison."""
    from migec.checkout import run as checkout_run
    from migec.sheet import preset

    out.mkdir(parents=True, exist_ok=True)
    sheet = out / "barcodes.txt"
    sheet.write_text(f"PBMC\t{preset('10x-v2')[0]}\n")

    co = out / "co"
    c = checkout_run(r1, sheet, co, r2, threads=threads)
    if c["total"] != CELLRANGER_READ_PAIRS:
        raise SystemExit(
            f"checkout read {c['total']:,} read pairs, not the {CELLRANGER_READ_PAIRS:,} Cell "
            f"Ranger reports. L001 alone is 3,155,166 -- concatenate BOTH lanes, in the same order "
            f"for R1 and R2, or the comparison is of two different libraries"
        )

    # Never: select the payload mate BY NAME. Paired checkout always writes `PBMC_R1.fq.gz` too,
    # holding one zero-length record per pair -- the min_payload gate is charged against the pair,
    # so a 26-of-26 trim on R1 is not a rejection. A glob picks up both, and refine on the empty
    # mate succeeds, reports clonality 1.0 and corrects on counts alone.
    mates = sorted(co.glob("PBMC_R2.fq.gz"))
    if len(mates) != 1:
        raise SystemExit(f"{len(mates)} files match {co}/PBMC_R2.fq.gz -- expected exactly one "
                         f"payload mate; refine on the barcode mate would report vacuous evidence")
    return c, mates[0]


def run_refine(mate: Path, out: Path, whitelist: Path | str, expect_cells: int,
               threads: int) -> dict:
    """Correct and call cells, once per whitelist setting."""
    from migec.refine import run as refine_run

    r = refine_run(mate, out, cell_whitelist=whitelist, expect_cells=expect_cells,
                   threads=threads)
    reads, total_reads, molecules = per_cell_reads(out / "PBMC.barcodes.tsv")
    return {
        "refine": r,
        "reads_per_cell": reads,
        "reads_total": total_reads,
        "molecules": molecules,
        "cells": called_cells(out / "PBMC.cells.tsv"),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", type=Path, required=True, help="R1, BOTH lanes concatenated")
    ap.add_argument("--r2", type=Path, required=True, help="R2, the same lanes in the same order")
    ap.add_argument("--cellranger-dir", type=Path, required=True,
                    help="directory holding 10x's published Cell Ranger outputs")
    ap.add_argument("--whitelist", type=Path, required=True, help="737K-august-2016.txt")
    ap.add_argument("--expect-cells", type=int, default=3000, help="OrdMag prior")
    ap.add_argument("--cellranger-version", default="5.0.0",
                    help="labels the Cell Ranger row; the published calls are 5.0.0")
    # The guard below exists to catch a wrong directory, not to pin a version -- so the expected
    # counts are settable, and having to set them is what makes pointing at another version
    # deliberate. 5.0.0 published 479 cells over 943 contigs; 10.1.0 run here gives 478 over 939.
    ap.add_argument("--expect-cellranger-cells", type=int, default=CELLRANGER_CELLS)
    ap.add_argument("--expect-cellranger-contigs", type=int, default=CELLRANGER_CONTIGS)
    ap.add_argument("--cellranger-seconds", type=float, default=None,
                    help="wall clock of a `cellranger vdj` run, from /usr/bin/time -v")
    ap.add_argument("--cellranger-rss-bytes", type=int, default=None,
                    help="peak RSS of that run, from /usr/bin/time -v (which reports kbytes)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--tsv", type=Path, default=None)
    args = ap.parse_args(argv)

    d = args.cellranger_dir.expanduser()
    contig_csv = next(iter(sorted(list(d.glob("*filtered_contig_annotations.csv")))), None)
    airr_tsv = next(iter(sorted(list(d.glob("*airr_rearrangement.tsv")))), None)
    metrics_csv = next(iter(sorted(list(d.glob("*metrics_summary.csv")))), None)
    for name, p in (("filtered_contig_annotations.csv", contig_csv),
                    ("airr_rearrangement.tsv", airr_tsv),
                    ("metrics_summary.csv", metrics_csv)):
        if p is None:
            raise SystemExit(f"no *_{name} in {d} -- see the header for the fetch command")

    check_airr_counts(airr_tsv, contig_csv)
    cr_cells, cr_contigs = read_cellranger_cells(contig_csv)
    if (len(cr_cells) != args.expect_cellranger_cells
            or cr_contigs != args.expect_cellranger_contigs):
        raise SystemExit(
            f"{len(cr_cells)} barcodes with is_cell true over {cr_contigs} contigs, not the "
            f"{args.expect_cellranger_cells} / {args.expect_cellranger_contigs} expected. Either "
            f"this is a different sample, or it is a different Cell Ranger version -- 5.0.0 "
            f"published 479/943 and 10.1.0 gives 478/939 on the same reads. Pass "
            f"--expect-cellranger-cells / --expect-cellranger-contigs to say which you meant"
        )
    metrics = read_metrics(metrics_csv)
    wl = read_whitelist(args.whitelist.expanduser())
    missing = sum(1 for b in cr_cells if b not in wl)
    if missing:
        raise SystemExit(f"{missing} of {len(cr_cells)} Cell Ranger cell barcodes are not on "
                         f"{args.whitelist} -- that is the wrong whitelist for this chemistry")

    args.out.mkdir(parents=True, exist_ok=True)
    rows = []

    # Cell Ranger's row. Only the three quantities it publishes on the whole library; its molecule
    # count is per contig and is deliberately absent (see the Never: in the header).
    rows.append({
        "tool": f"cellranger-{args.cellranger_version}",
        "whitelist": "737K-august-2016",
        "expect_cells": "",
        "reads_total": CELLRANGER_READ_PAIRS,
        "frac_valid_barcode": _pct(metrics["Valid Barcodes"]),
        "cells_called": len(cr_cells),
        "cells_shared": "",
        "cells_tool_only": "",
        "cells_other_only": "",
        "molecules": "",
        "frac_reads_in_cells": _pct(metrics["Fraction Reads in Cells"]),
        "wall_seconds": args.cellranger_seconds if args.cellranger_seconds is not None else "",
        "peak_rss_bytes": args.cellranger_rss_bytes if args.cellranger_rss_bytes is not None else "",
    })

    # Demultiplex ONCE -- the whitelist is a refine knob. Its clock is charged to BOTH migec rows,
    # because Cell Ranger's single number covers demultiplexing too.
    t0 = time.perf_counter()
    _, mate = run_checkout(args.r1, args.r2, args.out, args.threads)
    checkout_seconds = time.perf_counter() - t0

    # Then refine twice: without the whitelist so the raw barcodes are visible and the
    # read-weighted validity can be measured against the list, and with it so the cell set is the
    # corrected one.
    for label, whitelist in (("off", ""), ("on", str(args.whitelist.expanduser()))):
        t1 = time.perf_counter()
        m = run_refine(mate, args.out / f"wl_{label}", whitelist, args.expect_cells, args.threads)
        refine_seconds = time.perf_counter() - t1
        # Read-weighted, because Cell Ranger's 90.6% is a share of READS. refine's own whitelist
        # tally counts DISTINCT BARCODES (include/migec/whitelist.hpp:47-53) -- only
        # `reads_corrected` is read-weighted -- so the read weighting here is this script's.
        valid = sum(n for cell, n in m["reads_per_cell"].items() if cell in wl)
        shared = m["cells"] & cr_cells
        in_cells = sum(m["reads_per_cell"].get(c, 0) for c in m["cells"])
        rows.append({
            "tool": "migec",
            "whitelist": "737K-august-2016" if whitelist else "none",
            "expect_cells": args.expect_cells,
            "reads_total": m["reads_total"],
            "frac_valid_barcode": valid / max(m["reads_total"], 1),
            "cells_called": len(m["cells"]),
            "cells_shared": len(shared),
            "cells_tool_only": len(m["cells"] - cr_cells),
            "cells_other_only": len(cr_cells - m["cells"]),
            "molecules": m["molecules"],
            "frac_reads_in_cells": in_cells / max(m["reads_total"], 1),
            "wall_seconds": checkout_seconds + refine_seconds,
            "peak_rss_bytes": m["refine"].get("peak_rss_bytes", ""),
        })
        w = m["refine"]["whitelist"]
        print(f"# migec whitelist {label}: {w['exact']:,} barcodes exact, {w['corrected']:,} "
              f"snapped, {w['off_list']:,} off-list, background prior "
              f"{w['background_prior']:.2e}", file=sys.stderr)

    cols = ["tool", "whitelist", "expect_cells", "reads_total", "frac_valid_barcode",
            "cells_called", "cells_shared", "cells_tool_only", "cells_other_only",
            "molecules", "frac_reads_in_cells", "wall_seconds", "peak_rss_bytes"]
    fmt = {"frac_valid_barcode": ".4f", "frac_reads_in_cells": ".4f", "wall_seconds": ".1f"}
    lines = ["\t".join(cols)]
    for row in rows:
        lines.append("\t".join(
            "" if row[c] == "" else format(row[c], fmt.get(c, ".0f")) if isinstance(row[c], float)
            else str(row[c])
            for c in cols
        ))
    out = "\n".join(lines)
    print(out)
    if args.tsv:
        args.tsv.write_text(out + "\n")

    mig = [r for r in rows if r["tool"] == "migec"]
    if mig:
        best = mig[-1]
        print(f"# migec calls {best['cells_called']:,} cells against Cell Ranger's "
              f"{CELLRANGER_CELLS}, sharing {best['cells_shared']:,}; reads in cells "
              f"{best['frac_reads_in_cells']:.1%} against {_pct(metrics['Fraction Reads in Cells']):.1%}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
