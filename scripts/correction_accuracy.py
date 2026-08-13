#!/usr/bin/env python3
# 2026-08-13  Where does UMI error correction stop working, and why?
#
# `correct_umis` decides that barcode x is an error child of barcode y from three things: the
# distance (1 substitution), the estimated per-base error rate, and the COUNT RATIO -- a child
# should be much smaller than its parent, and `max_child_fraction` refuses a merge outright when
# it is not.
#
# That last piece is the whole game on a deeply sequenced amplicon and it is *absent* on a library
# sequenced at 1-3 reads per UMI, which is what bulk repertoire profiling and shallow 3' GEX both
# look like. A parent with 2 reads and a child with 1 is not an asymmetry; two singletons carry no
# count evidence whatsoever.
#
# So measure it rather than assume it. The simulator records `umi_true` and `umi_observed` per
# read, so an error child is exactly an observed barcode that is never its own true barcode, and
# recall/precision are countable against that.
#
#     python scripts/correction_accuracy.py

from __future__ import annotations

import argparse
import collections
import pathlib
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def evaluate(cov, molecules, umi_error, umi_len, out_dir, seed=5):
    """Simulate at one depth, correct, and score against the truth files."""
    from tests.synthetic._sim import SimConfig, simulate

    from migec import _core

    sim = simulate(
        SimConfig(adapter=ADAPTER, n_molecules=molecules, n_clones=100, coverage=cov,
                  coverage_cv=0.4, umi_len=umi_len, umi_error=umi_error, seed=seed),
        out_dir / f"sim{cov}",
    )
    parents = collections.defaultdict(collections.Counter)
    observed = collections.Counter()
    with open(sim["truth_reads"]) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            parents[f[4]][f[3]] += 1
            observed[f[4]] += 1

    # An error child is an observed barcode that is never its own true barcode. Anything else is a
    # real molecule, and merging one away destroys a molecule that was really there.
    children = {u for u, t in parents.items() if u not in t}
    real = [u for u in observed if u not in children]

    result = _core.correct_umis(list(observed.elements()))
    merged = {m["from"]: m["to"] for m in result["merges"]}
    correct = sum(1 for c, p in merged.items() if c in children and p in parents[c])
    return {
        "coverage": cov,
        "reads_per_umi": sum(observed.values()) / len(observed),
        "children": len(children),
        "merged": len(merged),
        "recall": correct / max(len(children), 1),
        "precision": correct / max(len(merged), 1),
        "molecules_kept": sum(1 for u in real if u not in merged) / max(len(real), 1),
        "epsilon": result["estimated_error"],
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--molecules", type=int, default=20_000)
    p.add_argument("--umi-error", type=float, default=3e-3)
    p.add_argument("--umi-length", type=int, default=12)
    p.add_argument("--coverage", type=float, nargs="+",
                   default=[1.3, 2.0, 3.0, 4.0, 6.0, 10.0, 25.0])
    p.add_argument("--out")
    a = p.parse_args(argv)

    out = pathlib.Path(a.out) if a.out else pathlib.Path(tempfile.mkdtemp())
    out.mkdir(parents=True, exist_ok=True)

    print(f"injected UMI error {a.umi_error:.1e} per base over {a.umi_length} nt, "
          f"{a.molecules:,} molecules\n")
    print(f"{'coverage':>9}{'reads/UMI':>10}{'children':>10}{'merged':>9}{'recall':>8}"
          f"{'precision':>11}{'molecules kept':>16}{'eps est':>10}{'eps/true':>10}")
    rows = []
    for cov in a.coverage:
        r = evaluate(cov, a.molecules, a.umi_error, a.umi_length, out)
        rows.append(r)
        print(f"{r['coverage']:>9.1f}{r['reads_per_umi']:>10.2f}{r['children']:>10,}"
              f"{r['merged']:>9,}{r['recall']:>8.3f}{r['precision']:>11.3f}"
              f"{r['molecules_kept']:>16.3f}{r['epsilon']:>10.2e}"
              f"{r['epsilon'] / a.umi_error:>10.2f}")

    with open(out / "correction_accuracy.tsv", "w") as fh:
        fh.write("coverage\treads_per_umi\tchildren\tmerged\trecall\tprecision\t"
                 "molecules_kept\tepsilon\n")
        for r in rows:
            fh.write(f"{r['coverage']}\t{r['reads_per_umi']:.4f}\t{r['children']}\t{r['merged']}\t"
                     f"{r['recall']:.4f}\t{r['precision']:.4f}\t{r['molecules_kept']:.4f}\t"
                     f"{r['epsilon']:.6e}\n")

    # The one number that decides whether the count-ratio rule is usable on a given library.
    usable = [r for r in rows if r["recall"] >= 0.8 and r["precision"] >= 0.9]
    if usable:
        print(f"\nrecall >= 0.8 and precision >= 0.9 from {min(r['reads_per_umi'] for r in usable):.1f} "
              f"reads/UMI upward.")
    print("Below that the count ratio carries no evidence -- a parent with 2 reads and a child")
    print("with 1 is not an asymmetry, and two singletons are not one either. What is still")
    print("available at 1 read is the barcode's own base QUALITY and the read's PAYLOAD: an error")
    print("child is a read of the parent's molecule, so its payload matches; an independent")
    print("molecule at distance 1 does not. Neither is used yet.")
    print(f"\nwrote {out / 'correction_accuracy.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
