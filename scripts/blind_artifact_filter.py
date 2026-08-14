#!/usr/bin/env python3
# 2026-08-14
# Tell a real low-frequency variant from a systematic artifact WITHOUT a panel of normals.
#
# The problem: a standard caller on UMI consensus reads emits systematic false positives. Measured
# on a 0%-certified reference arm, 91-94% of them were `-> G` -- the signature of 2-colour
# chemistry, where G is the base call for NO SIGNAL -- and eight positions recurred in every
# replicate. A consensus cannot remove them, because they are context-driven rather than random:
# most reads of the molecule read dark at the same base.
#
# The usual fix is a per-position background model from a panel of normals. You do not always have
# one. This is what you can do blind, and it works because migec knows something a caller does not:
# HOW MANY READS SUPPORTED EACH MOLECULE.
#
#   A real variant at VAF f  sits in ~f of molecules at full within-molecule support, and does not
#                            care how many reads built the molecule.
#   A context artifact       is carried disproportionately by molecules built from FEW reads --
#                            a singleton "consensus" IS one raw read, with no error correction at
#                            all, so it inherits the raw per-base error rate directly.
#
# So call the same sample at several `--min-reads` thresholds and watch each variant:
#
#   VAF stable as the threshold rises  -> real
#   VAF collapses, or the call vanishes -> it was carried by uncorrected small families
#
# Measured on certified material at 20 ng / 10x: every `-> G` artifact disappeared at
# `--min-reads 3`, and every certified variant survived with VAF stable to the third decimal.
#
# Usage:
#   migec assemble rf/S1.fq.gz -o as1/ --min-reads 1
#   migec assemble rf/S1.fq.gz -o as3/ --min-reads 3
#   migec assemble rf/S1.fq.gz -o as5/ --min-reads 5
#   ... align and call each, then:
#   python scripts/blind_artifact_filter.py --vcf 1=s1.mr1.vcf 3=s1.mr3.vcf 5=s1.mr5.vcf

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

# Below this the call is treated as absent rather than as a measured zero: a caller that stops
# reporting a site is telling you it lost significance, not that the frequency became 0.
ABSENT = None


def read_vcf(path: Path) -> dict[tuple[str, str, str, str], float]:
    """(chrom, pos, ref, alt) -> allele frequency, from a VCF or VCF.gz."""
    opener = gzip.open if path.suffix == ".gz" else open
    out: dict[tuple[str, str, str, str], float] = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            m = re.search(r"(?:^|;)AF=([0-9.eE+-]+)", f[7])
            if m:
                out[(f[0], f[1], f[3], f[4])] = float(m.group(1))
    return out


def classify(vafs: dict[int, float | None], drop: float, min_seen: int) -> tuple[str, float]:
    """Call it real or artifact from how the VAF moves with the min-reads threshold.

    `ratio` is the highest threshold's VAF over the lowest's. A variant that vanishes scores 0.
    Never: a call missing at a HIGHER threshold is evidence, but a call missing at the LOWEST one
    is not -- that is a caller gaining power, which is the opposite situation, so it is reported
    as inconclusive rather than being scored.
    """
    order = sorted(vafs)
    lo, hi = order[0], order[-1]
    if vafs[lo] is None:
        return "inconclusive", float("nan")
    seen = sum(1 for v in vafs.values() if v is not None)
    if seen < min_seen:
        return "artifact", 0.0
    if vafs[hi] is None:
        return "artifact", 0.0
    ratio = vafs[hi] / vafs[lo] if vafs[lo] else float("inf")
    return ("real" if ratio >= drop else "artifact"), ratio


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", nargs="+", required=True, metavar="MINREADS=PATH",
                    help="one per threshold, e.g. 1=a.vcf 3=b.vcf 5=c.vcf")
    ap.add_argument("--drop", type=float, default=0.5,
                    help="keep a variant if VAF at the highest threshold is at least this "
                         "fraction of the lowest (default 0.5)")
    ap.add_argument("--min-seen", type=int, default=2,
                    help="thresholds a variant must still be called at (default 2)")
    ap.add_argument("--out", type=Path, help="write the table here as TSV")
    a = ap.parse_args()

    tables: dict[int, dict] = {}
    for spec in a.vcf:
        if "=" not in spec:
            raise SystemExit(f"--vcf takes MINREADS=PATH, got {spec!r}")
        k, p = spec.split("=", 1)
        tables[int(k)] = read_vcf(Path(p))
    if len(tables) < 2:
        raise SystemExit("need at least two thresholds to see a trend")

    thresholds = sorted(tables)
    keys = sorted({k for t in tables.values() for k in t}, key=lambda k: (k[0], int(k[1])))

    rows = []
    for key in keys:
        vafs = {mr: tables[mr].get(key, ABSENT) for mr in thresholds}
        verdict, ratio = classify(vafs, a.drop, a.min_seen)
        rows.append({
            "chrom": key[0], "pos": key[1], "ref": key[2], "alt": key[3],
            **{f"vaf_mr{mr}": ("" if vafs[mr] is None else f"{vafs[mr]:.6f}") for mr in thresholds},
            "ratio": "" if ratio != ratio else f"{ratio:.2f}",   # nan-safe
            "verdict": verdict,
        })

    cols = list(rows[0]) if rows else []
    width = {c: max(len(c), 9) for c in cols}
    print("  ".join(c.ljust(width[c]) for c in cols))
    for r in rows:
        print("  ".join(str(r[c]).ljust(width[c]) for c in cols))

    real = [r for r in rows if r["verdict"] == "real"]
    art = [r for r in rows if r["verdict"] == "artifact"]
    print(f"\n{len(real)} real, {len(art)} artifact, "
          f"{len(rows) - len(real) - len(art)} inconclusive")
    if art:
        subs = [f"{r['ref']}>{r['alt']}" for r in art]
        tog = sum(1 for s in subs if s.endswith(">G"))
        print(f"artifact substitutions: {tog}/{len(subs)} are ->G"
              + ("  (2-colour dark-G signature)" if tog / len(subs) > 0.6 else ""))

    if a.out:
        a.out.write_text("\t".join(cols) + "\n"
                         + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows))
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
