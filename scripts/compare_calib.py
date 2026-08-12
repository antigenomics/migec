#!/usr/bin/env python3
# 2026-08-13
#
# UMI grouping accuracy against Calib (Orabi et al., Bioinformatics 2019; github.com/vpc-ccg/calib).
#
# Both tools answer the same question -- which reads came from the same original molecule -- so the
# comparison is a clustering comparison, scored against a known truth with the adjusted Rand index.
# Calib clusters on barcode AND read sequence with a locality-sensitive index; migec (today) groups
# on the barcode alone, which is the honest state of the pipeline until `assemble` lands and can
# split a group by its sequence. Expect Calib to win on libraries where distinct molecules collide
# on a barcode, and expect the gap to be the size of the collision rate.
#
# Usage:
#     python scripts/compare_calib.py --truth truth_reads.tsv \
#         --migec out/S1.fq.gz [--migec out/S2.fq.gz ...] \
#         --calib calib_out.cluster
#
# `truth_reads.tsv` is `read_id<TAB>molecule_id` with a header; tests/synthetic/_sim.py writes one
# directly. Any number of tools may be omitted -- the script reports whatever it was given.

from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path


def read_truth(path: Path) -> dict[str, str]:
    """read_id -> molecule_id. Extra columns are ignored, so the simulator's file works as is."""
    out: dict[str, str] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_read = header.index("read_id")
            i_mol = header.index("molecule_id")
        except ValueError:
            raise SystemExit(f"{path}: need 'read_id' and 'molecule_id' columns, got {header}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[i_read]] = f[i_mol]
    return out


def read_migec(paths: list[Path]) -> dict[str, str]:
    """read_id -> group, taken from the RX/BC tags checkout writes into the FASTQ header."""
    out: dict[str, str] = {}
    for path in paths:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4:
                    continue
                name, _, comment = line[1:].rstrip("\n").partition(" ")
                tags = dict(
                    (t.split(":", 2)[0], t.split(":", 2)[2])
                    for t in comment.split("\t")
                    if t.count(":") >= 2
                )
                if "RX" not in tags:
                    raise SystemExit(f"{path}: read {name} has no RX tag -- was this checkout out?")
                out[name] = f"{tags.get('BC', '')}:{tags['RX']}"
    return out


def read_calib(path: Path) -> dict[str, str]:
    """read_id -> cluster, from Calib's .cluster file.

    Nine tab-separated columns, one line per read: cluster_id, node_id, read_id, then name/seq/qual
    for each mate. Both mate names map to the same cluster; the forward one is used, with its
    trailing /1 stripped so it lines up with the truth file.
    """
    out: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            name = f[3].split()[0]
            for suffix in ("/1", "/2"):
                if name.endswith(suffix):
                    name = name[: -len(suffix)]
            out[name] = f[0]
    return out


def adjusted_rand(truth: dict[str, str], pred: dict[str, str]) -> dict[str, float]:
    """ARI plus the two numbers that say *how* a partition is wrong.

    ARI alone hides the direction of the error, and the two directions have opposite costs:
    over-splitting inflates the molecule count, over-merging destroys real variants. Both are
    reported, as the fraction of reads whose truth-group is split across predicted clusters and the
    fraction whose predicted cluster mixes truth-groups.
    """
    shared = [r for r in truth if r in pred]
    if not shared:
        raise SystemExit("no read names in common between the truth and the prediction")

    contingency: dict[tuple[str, str], int] = defaultdict(int)
    a: dict[str, int] = defaultdict(int)
    b: dict[str, int] = defaultdict(int)
    for r in shared:
        contingency[(truth[r], pred[r])] += 1
        a[truth[r]] += 1
        b[pred[r]] += 1

    n = len(shared)

    def c2(x: int) -> float:
        return x * (x - 1) / 2.0

    sum_ij = sum(c2(v) for v in contingency.values())
    sum_a = sum(c2(v) for v in a.values())
    sum_b = sum(c2(v) for v in b.values())
    expected = sum_a * sum_b / c2(n) if n > 1 else 0.0
    maximum = 0.5 * (sum_a + sum_b)
    ari = (sum_ij - expected) / (maximum - expected) if maximum != expected else 1.0

    # Reads sitting in a truth group that ended up spread over more than one predicted cluster...
    per_truth_clusters: dict[str, set[str]] = defaultdict(set)
    per_pred_truths: dict[str, set[str]] = defaultdict(set)
    for t, p in contingency:
        per_truth_clusters[t].add(p)
        per_pred_truths[p].add(t)
    split = sum(a[t] for t in per_truth_clusters if len(per_truth_clusters[t]) > 1)
    merged = sum(b[p] for p in per_pred_truths if len(per_pred_truths[p]) > 1)

    return {
        "reads": float(n),
        "true_molecules": float(len(a)),
        "predicted_clusters": float(len(b)),
        "ari": ari,
        "reads_in_split_molecules": split / n,
        "reads_in_merged_clusters": merged / n,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", type=Path, required=True, help="read_id/molecule_id TSV")
    ap.add_argument("--migec", type=Path, action="append", default=[],
                    help="checkout FASTQ (repeatable, one per sample)")
    ap.add_argument("--calib", type=Path, help="Calib .cluster file")
    ap.add_argument("--partition", type=Path, action="append", default=[],
                    help="any other tool, as a read_id/cluster_id TSV named tool=path")
    args = ap.parse_args(argv)

    truth = read_truth(args.truth)
    results: dict[str, dict[str, float]] = {}
    if args.migec:
        results["migec"] = adjusted_rand(truth, read_migec(args.migec))
    if args.calib:
        results["calib"] = adjusted_rand(truth, read_calib(args.calib))
    for spec in args.partition:
        name, _, path = str(spec).partition("=")
        pred = {}
        with open(path or name) as fh:
            for line in fh:
                r, _, c = line.rstrip("\n").partition("\t")
                pred[r] = c
        results[name] = adjusted_rand(truth, pred)

    if not results:
        raise SystemExit("nothing to compare: pass at least one of --migec/--calib/--partition")

    cols = ["reads", "true_molecules", "predicted_clusters", "ari",
            "reads_in_split_molecules", "reads_in_merged_clusters"]
    print("tool\t" + "\t".join(cols))
    for tool, r in results.items():
        cells = [f"{r[c]:.0f}" if r[c] >= 1000 or c.endswith(("reads", "molecules", "clusters"))
                 else f"{r[c]:.4f}" for c in cols]
        print(tool + "\t" + "\t".join(cells))
    return 0


if __name__ == "__main__":
    sys.exit(main())
