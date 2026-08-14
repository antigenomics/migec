#!/usr/bin/env python3
# 2026-08-14
# Per-TARGET molecule counts and variant calls for the ctDNA titration, against certified truth.
#
# This is the reference-based answer to what `ctdna_titration.py` could only average. That script
# divides a library total by an amplicon count inferred from consensus prefixes; this one takes
# the molecules actually aligned to each target of a panel inferred from coverage, so
# "molecules at a site" is a count rather than a mean.
#
# Inputs are produced by `integrations/../job/per_site.sbatch` on a cluster:
#   molecules_per_target.tsv   run, chrom, start, end, gene, molecules
#   variants.tsv               run, chrom, pos, ref, alt, dp, af
#   design.tsv                 run_accession, sample_alias
#
# Usage:
#   python scripts/ctdna_persite.py --molecules m.tsv --variants v.tsv --design d.tsv --out out/

from __future__ import annotations

import argparse
import collections
import math
import statistics as st
from pathlib import Path

# The certified mutant allele frequencies of the reference material, from the sample alias.
ARMS = {"WT": 0.0, "0.125": 0.00125, "0.25": 0.0025, "1": 0.01}
# The undiluted material's frequency is not certified, so it is carried as unknown.
UNDILUTED = "cell_line"


def parse_alias(alias: str) -> dict:
    """`<input>ng_<arm>_<depth>x_rep_<n>` -> the design fields, or empties."""
    out = {"input_ng": None, "arm": None, "vaf": None, "depth": None, "rep": None,
           "diluted": None}
    parts = alias.split("_")
    if len(parts) < 4 or not parts[0].endswith("ng"):
        return out
    out["input_ng"] = int(parts[0][:-2])
    arm = parts[1] if parts[1] != UNDILUTED else UNDILUTED
    if parts[1] == "cell" and len(parts) > 2 and parts[2] == "line":
        arm = UNDILUTED
        rest = parts[3:]
    else:
        rest = parts[2:]
    out["arm"] = arm
    out["vaf"] = ARMS.get(arm)
    # Never: the undiluted arm is NOT a 0% control and NOT a certified frequency. It is the raw
    # material, and treating it as either would put a fabricated number into the comparison.
    out["diluted"] = arm != UNDILUTED
    if rest and rest[0].endswith("x"):
        out["depth"] = float(rest[0][:-1])
    if "rep" in parts:
        out["rep"] = int(parts[parts.index("rep") + 1])
    return out


def read_tsv(path: Path) -> list[dict]:
    lines = path.read_text().strip().split("\n")
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:] if ln]


def poisson_at_least(lam: float, k: int) -> float:
    if lam <= 0:
        return 0.0
    return 1 - sum(math.exp(-lam) * lam**i / math.factorial(i) for i in range(k))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--molecules", type=Path, required=True)
    ap.add_argument("--variants", type=Path, required=True)
    ap.add_argument("--design", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--hotspot", default="3:179234297",
                    help="the certified variant's locus (default PIK3CA H1047R)")
    ap.add_argument("--min-support", type=int, default=3)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    design = {r["run_accession"]: r["sample_alias"] for r in read_tsv(a.design)}
    meta = {run: parse_alias(alias) for run, alias in design.items()}

    mols = collections.defaultdict(dict)
    for r in read_tsv(a.molecules):
        mols[r["run"]][f"{r['chrom']}:{r['start']}-{r['end']}"] = (int(r["molecules"]), r["gene"])

    calls = collections.defaultdict(list)
    for r in read_tsv(a.variants):
        calls[r["run"]].append((f"{r['chrom']}:{r['pos']}", f"{r['ref']}>{r['alt']}",
                                float(r["af"]), int(r["dp"])))

    runs = [r for r in sorted(mols) if meta.get(r, {}).get("input_ng")]

    # --- per-target molecule counts, and what the weakest target costs -------------------------
    rows = []
    for run in runs:
        m = meta[run]
        per = {k: v[0] for k, v in mols[run].items()}
        if not per:
            continue
        vals = sorted(per.values())
        rows.append({
            "run": run, "arm": m["arm"], "input_ng": m["input_ng"], "depth": m["depth"],
            "targets": len(per), "total": sum(vals), "mean": round(st.mean(vals)),
            "weakest": vals[0], "strongest": vals[-1],
            "weakest_over_mean": round(vals[0] / st.mean(vals), 3),
            "calls": len(calls.get(run, [])),
        })
    hdr = list(rows[0])
    (a.out / "per_target.tsv").write_text(
        "\t".join(hdr) + "\n" + "".join("\t".join(str(r[c]) for c in hdr) + "\n" for r in rows))

    ratios = [r["weakest_over_mean"] for r in rows]
    print(f"{len(rows)} runs, {rows[0]['targets']} targets")
    print(f"weakest target holds {min(ratios):.2f}-{max(ratios):.2f}x the panel mean "
          f"(median {st.median(ratios):.2f})")

    # --- specificity: what the 0%-certified arm reports ------------------------------------------
    print("\ncall burden by arm (the specificity view):")
    burden = collections.defaultdict(list)
    for r in rows:
        burden[r["arm"]].append(r["calls"])
    for arm in sorted(burden, key=lambda x: (x == UNDILUTED, x)):
        v = burden[arm]
        print(f"  {arm:<10} n={len(v):<3} {st.mean(v):>5.1f} calls/sample  "
              f"range {min(v)}-{max(v)}")

    # --- the test that separates the two explanations --------------------------------------------
    # Does artifact burden track MOLECULE COUNT (statistical power) or SAMPLE PREPARATION
    # (dilution handling)? They are confounded in the raw design, so compare at matched molecules.
    print("\ncall burden against molecules, split by preparation:")
    print(f"  {'preparation':<14}{'molecules/site':>16}{'calls/sample':>14}{'n':>5}")
    for diluted in (True, False):
        sel = [r for r in rows if meta[r["run"]]["diluted"] is diluted]
        if not sel:
            continue
        lo = [r for r in sel if r["mean"] < st.median([x["mean"] for x in sel])]
        hi = [r for r in sel if r["mean"] >= st.median([x["mean"] for x in sel])]
        for label, grp in (("low", lo), ("high", hi)):
            if grp:
                print(f"  {'diluted' if diluted else 'undiluted':<9}{label:<5}"
                      f"{st.mean([r['mean'] for r in grp]):>16,.0f}"
                      f"{st.mean([r['calls'] for r in grp]):>14.1f}{len(grp):>5}")

    # --- the hotspot, arm by arm -------------------------------------------------------------------
    print(f"\ncertified locus {a.hotspot}, measured against truth:")
    print(f"  {'arm':<12}{'certified':>11}{'n':>4}{'detected':>10}{'mean VAF':>11}")
    by_arm = collections.defaultdict(list)
    for run in runs:
        m = meta[run]
        af = next((c[2] for c in calls.get(run, []) if c[0] == a.hotspot), 0.0)
        by_arm[m["arm"]].append(af)
    for arm in sorted(by_arm, key=lambda x: (x == UNDILUTED, x)):
        v = by_arm[arm]
        cert = ARMS.get(arm)
        det = sum(1 for x in v if x > 0)
        print(f"  {arm:<12}{(f'{cert:.5f}' if cert is not None else 'n/a'):>11}{len(v):>4}"
              f"{f'{det}/{len(v)}':>10}{st.mean(v):>11.4f}")

    # --- the substitution spectrum of the true negative ---------------------------------------------
    neg = [c for run in runs if meta[run]["arm"] == "WT" for c in calls.get(run, [])]
    if neg:
        spec = collections.Counter(c[1] for c in neg)
        tog = sum(v for k, v in spec.items() if k.endswith(">G"))
        print(f"\ntrue-negative substitution spectrum ({len(neg)} calls): "
              f"{tog}/{len(neg)} = {tog/len(neg):.0%} are ->G")
        for k, v in spec.most_common():
            print(f"  {k}  {v}")
        rec = collections.Counter((c[0], c[1]) for c in neg)
        n_wt = sum(1 for r in runs if meta[r]["arm"] == "WT")
        systematic = [(k, v) for k, v in rec.items() if v >= max(2, n_wt // 2)]
        print(f"  positions recurring in >= half of {n_wt} WT runs: {len(systematic)}")

    print(f"\nwrote {a.out / 'per_target.tsv'}")


if __name__ == "__main__":
    main()
