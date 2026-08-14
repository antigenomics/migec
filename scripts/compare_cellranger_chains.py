#!/usr/bin/env python3
# 2026-08-14
#
# The per-cell chain axis of the Cell Ranger comparison: given the same droplet library, do migec
# plus arda put the same receptor chain in the same cell as Cell Ranger's assembler?
#
# The two pipelines get there differently and that is the point. Cell Ranger assembles a per-CELL
# contig from every read of a barcode, then annotates it. migec assembles a per-MOLECULE consensus
# (`assemble --contig`, one component per overlap group), annotates each molecule with arda, and
# lets the cell's chain call be a vote over its molecules. So this scores the same question --
# which chain is in which cell -- reached with and without a per-cell assembler.
#
# Never: `restrict to molecules deep enough to carry a junction, and report the denominator.`
# `molecules_scored` and `molecules_total` are columns so it is never implicit.
#
# Note: **depth does not buy coverage here, and the arithmetic that says it does is wrong.** The
# tempting model is that each read is placed uniformly over the ~508 nt amplicon, so a 90 nt read
# spans the median 42 nt junction with p = 0.114 and 30 reads give 1 - 0.886^30 = 0.975. Measured,
# that is false: reads of one (CB, UMI) are co-terminal in this chemistry (X1: 92% of groups), so a
# molecule is a PILE at one position, not a tiling of the amplicon, and its consensus covers one
# window however deep it is. At `--min-reads 30` the mean consensus is 204 nt, not 508, and
# **7,855 of 47,584 molecules (16.5%) carry a cell, a locus and a junction** -- close to the
# single-window 0.32 the geometry predicts and nowhere near 0.975. The depth cut is still right,
# because a deep pile gives a clean consensus; it just does not extend it.
#
# Never: `recall is the primary metric, not concordance.` Concordance is computed over the chains
# BOTH tools called, which is a self-selecting denominator that reads high whatever happens. A
# missed chain is the unrecoverable error, so `chain_recall = chains_shared / chains_cellranger`
# leads and concordance follows it.
#
# Never: `arda runs with --exact, not the amplicon preset.` `amplicon` turns on `two_pass`,
# `fast_segments` and `segment_only_v`, which assume a primer-anchored read spanning V into J.
# These are 90 nt tiles of a fragmented amplicon -- the documented loss case -- so the preset would
# confound the comparison with a preset choice. `arda rnaseq --exact` is the zero preset.
#
# Never: `Cell Ranger's AIRR count columns are the inverse of arda's.` `consensus_count` is READS
# and `duplicate_count` is UMIs there; arda spells them the other way round. Nothing here joins on
# those names -- the chain identity is (cell, locus, junction_aa).
#
# Usage:
#     python scripts/compare_cellranger_chains.py \
#         --consensus /tmp/cr/asm30/PBMC.fq.gz \
#         --cellranger-dir ~/data/10x/sc5p_v2_hs_PBMC_1k_t_cellranger \
#         --out /tmp/cr/chains --min-reads 30 --tsv assets/cellranger_chains.tsv
#
# The consensus comes from:
#     migec assemble <refined.fq.gz> -o asm/ --contig --min-reads 30

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from compare_cellranger import read_cellranger_cells  # noqa: E402


def cellranger_chains(contig_csv: Path) -> dict[tuple[str, str], set[str]]:
    """(cell, locus) -> the set of CDR3 amino-acid strings Cell Ranger called there.

    `is_cell` filtered, and the `-1` gem-group suffix stripped so the key matches migec's bare
    16 nt barcode. A contig with no CDR3 call contributes the cell but no sequence, which is why
    the value is a set and not a string.
    """
    out: dict[tuple[str, str], set[str]] = defaultdict(set)
    with open(contig_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["is_cell"].strip().lower() != "true":
                continue
            cdr3 = (row.get("cdr3") or "").strip()
            if cdr3 and cdr3 != "None":
                out[(row["barcode"].split("-", 1)[0], row["chain"])].add(cdr3)
    return dict(out)


def run_arda(consensus: Path, out_dir: Path, arda: str, threads: int) -> Path:
    """`arda rnaseq --exact --cell-from migec`, returning the per-read AIRR TSV."""
    out_dir.mkdir(parents=True, exist_ok=True)
    argv = [arda, "rnaseq", "--exact", "--cell-from", "migec",
            "--r1", str(consensus), "--out-dir", str(out_dir), "--out-prefix", "cells",
            "--threads", str(threads)]
    p = subprocess.run(argv, capture_output=True, text=True)
    if p.returncode:
        raise SystemExit(f"arda failed ({p.returncode}):\n{p.stderr[-2000:]}")
    airr = out_dir / "cells.airr.tsv"
    if not airr.exists():
        raise SystemExit(f"arda wrote no {airr}; it reported:\n{p.stderr[-2000:]}")
    return airr


def migec_chains(airr: Path, min_molecules: int) -> tuple[dict, dict]:
    """(cell, locus) -> the junction_aa a vote over that cell's molecules picks, and the votes.

    A cell's chain call is the modal junction over its molecules. Never the first molecule: the
    consensus file is written in barcode order, so "first" is a property of the barcode.
    """
    votes: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    with open(airr, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "cell_id" not in (reader.fieldnames or []):
            raise SystemExit(
                f"{airr} has no `cell_id` column -- arda ran without `--cell-from`, so every "
                f"molecule is cellless and every recall below would be 0.0 with nothing saying why"
            )
        for row in reader:
            cell, locus = row.get("cell_id", ""), row.get("locus", "")
            junction = (row.get("junction_aa") or "").strip()
            if not cell or not locus or not junction:
                continue
            votes[(cell, locus)][junction] += 1
    called = {k: max(v.items(), key=lambda kv: (kv[1], kv[0]))[0]
              for k, v in votes.items() if sum(v.values()) >= min_molecules}
    return called, {k: sum(v.values()) for k, v in votes.items()}


def count_consensus(path: Path) -> int:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return sum(1 for i, _ in enumerate(fh) if i % 4 == 0)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--consensus", type=Path, required=True,
                    help="migec assemble --contig output, already depth-restricted")
    ap.add_argument("--cellranger-dir", type=Path, required=True)
    ap.add_argument("--arda", default=shutil.which("arda") or "arda")
    ap.add_argument("--cellranger-version", default="5.0.0",
                    help="labels the rows; two runs of this script concatenate into one table, so "
                         "the version has to be a COLUMN and not something you remember")
    ap.add_argument("--min-reads", type=int, required=True,
                    help="the depth cut already applied to --consensus; reported, not applied")
    ap.add_argument("--min-molecules", type=int, default=1,
                    help="molecules a (cell, locus) needs before its chain is called")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--tsv", type=Path, default=None)
    args = ap.parse_args(argv)

    d = args.cellranger_dir.expanduser()
    contig_csv = next(iter(sorted(list(d.glob("*filtered_contig_annotations.csv")))), None)
    if contig_csv is None:
        raise SystemExit(f"no *_filtered_contig_annotations.csv in {d}")
    cr_cells, _ = read_cellranger_cells(contig_csv)
    cr = cellranger_chains(contig_csv)

    if not args.consensus.exists():
        raise SystemExit(f"{args.consensus} does not exist -- run `migec assemble --contig "
                         f"--min-reads {args.min_reads}` first")
    molecules_total = count_consensus(args.consensus)

    airr = run_arda(args.consensus, args.out, args.arda, args.threads)
    mig, mig_votes = migec_chains(airr, args.min_molecules)
    molecules_scored = sum(mig_votes.values())

    # Score only where BOTH tools could have called: Cell Ranger's cells. A migec cell Cell Ranger
    # never called is not a miss by either tool, it is the population difference the other table
    # already reports.
    loci = sorted({locus for (_, locus) in cr})
    rows = []
    for locus in loci:
        cr_here = {c: j for (c, ln), j in cr.items() if ln == locus and c in cr_cells}
        mig_here = {c: j for (c, ln), j in mig.items() if ln == locus and c in cr_cells}
        shared = set(cr_here) & set(mig_here)
        agree = sum(1 for c in shared if mig_here[c] in cr_here[c])
        rows.append({
            "cellranger_version": args.cellranger_version,
            "locus": locus,
            "min_reads": args.min_reads,
            "min_molecules": args.min_molecules,
            "cells_scored": len(cr_cells),
            "chains_cellranger": len(cr_here),
            "chains_migec": len(mig_here),
            "chains_shared": len(shared),
            "chain_recall": len(shared) / max(len(cr_here), 1),
            "junction_aa_concordance": agree / max(len(shared), 1),
            "molecules_scored": molecules_scored,
            "molecules_total": molecules_total,
        })

    cols = ["cellranger_version", "locus", "min_reads", "min_molecules", "cells_scored", "chains_cellranger",
            "chains_migec", "chains_shared", "chain_recall", "junction_aa_concordance",
            "molecules_scored", "molecules_total"]
    fmt = {"chain_recall": ".4f", "junction_aa_concordance": ".4f"}
    lines = ["\t".join(cols)]
    for row in rows:
        lines.append("\t".join(
            format(row[c], fmt.get(c, ".0f")) if isinstance(row[c], float) else str(row[c])
            for c in cols))
    out = "\n".join(lines)
    print(out)
    if args.tsv:
        args.tsv.write_text(out + "\n")

    print(f"# {molecules_total:,} consensuses at >= {args.min_reads} reads; "
          f"{molecules_scored:,} carried a cell, a locus and a junction", file=sys.stderr)
    for row in rows:
        print(f"# {row['locus']}: recall {row['chain_recall']:.1%} "
              f"({row['chains_shared']}/{row['chains_cellranger']}), "
              f"junction agreement {row['junction_aa_concordance']:.1%}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
