#!/usr/bin/env python3
# 2026-08-15
#
# The contig axis of the Cell Ranger comparison: with no germline reference in the assembly at
# all, does migec plus arda reconstruct the same per-cell V(D)J contigs Cell Ranger's assembler
# does?
#
# The two get there differently. Cell Ranger assembles a per-CELL contig from every read of a
# barcode and annotates it against a germline reference. migec collapses each MOLECULE to a
# consensus, and `arda cells` overlap-assembles a cell's molecules into contigs -- k-mer seeds,
# verified overlaps, union-find layout, weighted column consensus. No reference is read until the
# contigs are annotated, which is after the assembly is finished.
#
# Never: `the premise is measured before the method is`. Reads of one (CB, UMI) are co-terminal in
# this chemistry, so a molecule's consensus is a pile at one position and covers one window of the
# transcript however deep it is -- at `--min-reads 30` the mean consensus is 204 nt against a
# ~508 nt amplicon. What makes the contig recoverable is that DIFFERENT molecules start at
# different positions, so a cell's molecules tile the transcript. The ceiling column of this table
# is that premise, scored directly: the share of each Cell Ranger contig's 25-mers present
# anywhere in its own cell's raw molecules, before any assembly. It is 0.999.
#
# Never: `each ingredient is measured by removing it, in this table, on this data`. The rows are
# the ablations, and CDR3-exact reads 0.9894 with everything on, 0.9820 without the phasing and
# 0.9618 without the adapter trim; chain recall 0.9777 / 0.9714 / 0.9491. Note the adapter row's
# contig N50 is the HIGHEST of the three (580 against 536) and its accuracy is the lowest -- a
# longer contig built across an adapter is a longer wrong contig, so N50 is a description and
# never a score. Depth weighting is the third ingredient and has no flag to turn off; it is
# measured in `arda.singlecell`'s own module table, at 0.967 against 0.988.
#
# Never: `a doublet is two chains of the same LOCUS`. One TRA and one TRB in a cell is a normal
# paired T cell. Two TRB is the doublet signal, and until the phasing landed it was invisible by
# construction: the two chains share their constant region, so the layout puts them in one
# component and the column consensus averages their junctions into a sequence that is neither.
#
# Never: `Cell Ranger is not re-run here`. The comparator is 10x's published output (5.0.0) and
# the local 10.1.0 run; both are inputs. See docs/single_cell.rst for the cost axis.
#
# Usage:
#     python scripts/compare_cellranger_contigs.py \
#         --consensus /tmp/cr/asm_all/PBMC.consensus.fq.gz \
#         --cellranger-dir ~/data/10x/sc5p_v2_hs_PBMC_1k_t_cellranger \
#         --arda ~/vcs/code/arda/.venv/bin/arda \
#         --out /tmp/cr/contigs --tsv assets/cellranger_contigs.tsv
#
# The consensus comes from, with EVERY molecule kept -- a shallow molecule is one more window of
# the tiling and `--min-reads` throws the tiling away:
#     migec assemble <refined.fq.gz> -o asm_all/ --contig --min-reads 1

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

K = 25
_COMP = str.maketrans("ACGTN", "TGCAN")


def rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def read_fasta(path: Path) -> dict[str, list[str]]:
    """cell -> its sequences, for a FASTA named either `<cell>_contig_<n>` or `<cell>-1_contig_n`."""
    out: dict[str, list[str]] = defaultdict(list)
    name, buf = None, []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name].append("".join(buf))
                name = line[1:].split()[0].split("_contig")[0].split("-")[0]
                buf = []
            elif line.strip():
                buf.append(line.strip().upper())
    if name:
        out[name].append("".join(buf))
    return out


def cellranger_contigs(fasta: Path, annotations: Path) -> tuple[dict, dict]:
    """(cell -> contig sequences, cell -> [(chain, cdr3_nt, cdr3_aa)]), `is_cell` filtered."""
    contigs = read_fasta(fasta)
    calls: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    with open(annotations, newline="") as fh:
        for row in csv.DictReader(fh):
            if str(row.get("is_cell", "true")).strip().lower() != "true":
                continue
            nt = (row.get("cdr3_nt") or "").strip()
            aa = (row.get("cdr3") or "").strip()
            if nt and nt != "None":
                calls[row["barcode"].split("-")[0]].append((row["chain"], nt, aa))
    return contigs, calls


def read_migec_molecules(consensus: Path, cells: set[str]) -> dict[str, list[str]]:
    """cell -> its molecule consensus sequences, restricted to the called cells."""
    out: dict[str, list[str]] = defaultdict(list)
    with gzip.open(consensus, "rt") as fh:
        for header, seq, _plus, _qual in zip(*[iter(fh)] * 4):
            if header[:1] != "@":
                continue
            cell = header[1:].split(None, 1)[0].split(".")[1]
            if cell in cells:
                out[cell].append(seq.strip().upper())
    return out


def kmer_coverage(target: str, pool: set[str]) -> float:
    """Share of `target`'s k-mers present in `pool`. A single mismatch costs 25 of them, so this
    is a strict measure of sequence agreement and not merely of overlap."""
    kmers = [target[i:i + K] for i in range(len(target) - K + 1)]
    if not kmers:
        return 0.0
    return sum(1 for x in kmers if x in pool) / len(kmers)


def kmers_of(sequences: list[str]) -> set[str]:
    pool: set[str] = set()
    for seq in sequences:
        for orientation in (seq, rc(seq)):
            for i in range(len(orientation) - K + 1):
                pool.add(orientation[i:i + K])
    return pool


def score(cellranger: tuple[dict, dict], ours: dict[str, list[str]]) -> dict:
    """k-mer coverage of each Cell Ranger contig, and exact CDR3 containment, per chain."""
    contigs, calls = cellranger
    coverage: list[float] = []
    cdr3_hit: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for cell, targets in contigs.items():
        pool = kmers_of(ours.get(cell, []))
        joined = [o for seq in ours.get(cell, []) for o in (seq, rc(seq))]
        for target in targets:
            coverage.append(kmer_coverage(target, pool))
        for chain, nt, _aa in calls.get(cell, []):
            found = any(nt in o for o in joined)
            cdr3_hit[chain][0] += found
            cdr3_hit[chain][1] += 1
    total = [sum(v[0] for v in cdr3_hit.values()), sum(v[1] for v in cdr3_hit.values())]
    return {
        "contigs_scored": len(coverage),
        "kmer_coverage_mean": sum(coverage) / len(coverage) if coverage else float("nan"),
        "contigs_at_90pct": sum(1 for x in coverage if x >= 0.9),
        "cdr3_exact": total[0],
        "cdr3_total": total[1],
        "cdr3_exact_rate": total[0] / total[1] if total[1] else float("nan"),
        "per_chain": {k: v for k, v in sorted(cdr3_hit.items())},
    }


def run_arda(arda: str, consensus: Path, out: Path, cells_file: Path, reference: Path,
             extra: list[str], threads: int) -> dict:
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [arda, "cells", str(consensus), "-p", str(out), "--cells", str(cells_file),
           "--reference", str(reference), "--threads", str(threads), *extra]
    subprocess.run(cmd, check=True)
    return json.loads(Path(str(out) + ".arda.json").read_text())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--consensus", type=Path, required=True,
                    help="migec assemble --contig --min-reads 1 output")
    ap.add_argument("--cellranger-dir", type=Path, required=True)
    ap.add_argument("--arda", default=shutil.which("arda") or "arda")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tsv", type=Path)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--cellranger-version", default="5.0.0")
    args = ap.parse_args()

    directory = args.cellranger_dir.expanduser()
    fasta = next(iter(sorted(directory.glob("*filtered_contig.fasta"))), None)
    annotations = next(iter(sorted(directory.glob("*filtered_contig_annotations.csv"))), None)
    if fasta is None or annotations is None:
        raise SystemExit(
            f"{directory} holds no *filtered_contig.fasta and *filtered_contig_annotations.csv "
            f"pair -- the contig axis needs the SEQUENCES, not only the annotation table")

    contigs, calls = cellranger_contigs(fasta, annotations)
    if not contigs:
        raise SystemExit(f"{fasta} yielded no contigs; the barcode naming is not what is assumed")
    print(f"cellranger {args.cellranger_version}: {len(contigs)} cells, "
          f"{sum(len(v) for v in contigs.values())} contigs, "
          f"{sum(len(v) for v in calls.values())} CDR3 calls", file=sys.stderr)

    args.out.mkdir(parents=True, exist_ok=True)
    cells_file = args.out / "cells.txt"
    cells_file.write_text("\n".join(sorted(contigs)) + "\n")

    molecules = read_migec_molecules(args.consensus, set(contigs))
    if not molecules:
        raise SystemExit(
            f"none of Cell Ranger's {len(contigs)} barcodes appears in {args.consensus}. "
            f"The consensus names must be migec's <sample>.<cell>.<umi>, and the run must not "
            f"have been restricted by --min-reads")
    print(f"migec: {sum(len(v) for v in molecules.values())} molecules over "
          f"{len(molecules)} of those cells", file=sys.stderr)

    rows: list[dict] = []

    # The premise, scored before any assembly runs: is the contig already in the molecules?
    ceiling = score((contigs, calls), molecules)
    rows.append({"method": "migec molecules (no assembly)", "variant": "ceiling",
                 "cells": len(molecules), **{k: v for k, v in ceiling.items()
                                             if k != "per_chain"},
                 "contig_n50": 0, "chain_recall": float("nan"),
                 "doublet_candidates": 0, "wall_note": "k-mer ceiling, not a method"})

    variants = [("arda cells", "default", []),
                ("arda cells", "no phasing", ["--no-split"]),
                ("arda cells", "no adapter trim", ["--no-trim-adapter"])]
    for method, variant, extra in variants:
        prefix = args.out / ("arda_" + variant.replace(" ", "_"))
        report = run_arda(args.arda, args.consensus, prefix, cells_file, annotations,
                          extra, args.threads)
        assembled = read_fasta(Path(str(prefix) + ".contigs.fasta"))
        scored = score((contigs, calls), assembled)
        rows.append({
            "method": method, "variant": variant, "cells": report["cells"],
            **{k: v for k, v in scored.items() if k != "per_chain"},
            "contig_n50": report["contig_n50"],
            "chain_recall": report.get("reference_recall", float("nan")),
            "doublet_candidates": report["cells_doublet_candidate"],
            "wall_note": f"{report['contigs']} contigs from {report['molecules_in']} molecules",
        })
        if variant == "default":
            for chain, (hit, total) in sorted(scored["per_chain"].items()):
                print(f"  {chain}: CDR3 exact {hit}/{total} ({hit / total:.4f})", file=sys.stderr)

    columns = ("method", "variant", "cells", "contigs_scored", "kmer_coverage_mean",
               "contigs_at_90pct", "cdr3_exact", "cdr3_total", "cdr3_exact_rate",
               "contig_n50", "chain_recall", "doublet_candidates", "wall_note")
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(
            f"{row[c]:.4f}" if isinstance(row[c], float) else str(row[c]) for c in columns))
    table = "\n".join(lines) + "\n"
    print(table)
    if args.tsv:
        args.tsv.write_text(table)
        print(f"wrote {args.tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
