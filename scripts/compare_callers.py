#!/usr/bin/env python3
# 2026-08-14
#
# Caller against caller, on the SAME consensus BAM, over the certified ctDNA arms.
#
# `docs/detection.rst` reports that LoFreq on migec consensus returns 9-11 calls per sample on the
# 0%-certified arm, 94% of them `-> G`, with eight positions recurring across all three replicates.
# That reads as a 2-colour dark-G artifact in the reads -- but one caller cannot establish it. If a
# second caller on the identical alignment reports the same positions, the artifact is in the
# input; if it does not, it was the first caller's threshold. `scripts/ctdna_callers.sbatch` runs
# both and this script scores the result.
#
# The 0% arm is the number that matters. Its variant frequency is zero by construction, so every
# call there is that pipeline's own false-positive rate on real chemistry rather than on a
# simulation.
#
# Usage:
#     python scripts/compare_callers.py --variants callers_variants.tsv \
#         --runtime callers_runtime.tsv --out assets/ctdna_callers.tsv

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ctdna_persite import ARMS, read_tsv  # noqa: E402

# The certified hotspot: PIK3CA H1047R, the one locus whose frequency the vendor certifies.
HOTSPOT = ("3", 179234297)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", type=Path, required=True)
    ap.add_argument("--runtime", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--hotspot", default=f"{HOTSPOT[0]}:{HOTSPOT[1]}")
    a = ap.parse_args(argv)

    chrom, pos = a.hotspot.split(":")
    calls = [r for r in read_tsv(a.variants) if r.get("filtered", "PASS") in ("PASS", ".")]
    if not calls:
        raise SystemExit(f"no PASS calls in {a.variants}")

    runs = collections.defaultdict(set)
    burden = collections.defaultdict(list)
    hotspot = collections.defaultdict(list)
    spectrum = collections.Counter()

    per_run = collections.defaultdict(list)
    for c in calls:
        key = (c["caller"], int(c["mapq_floor"]), int(c["min_reads"]), c["arm"])
        per_run[(key, c["run"])].append(c)
        runs[key].add(c["run"])

    for (key, run), cs in per_run.items():
        burden[key].append(len(cs))
        hit = [c for c in cs if c["chrom"] == chrom and c["pos"] == pos]
        hotspot[key].append(float(hit[0]["af"]) if hit else 0.0)
        if key[3] == "WT":
            for c in cs:
                spectrum[(c["ref"], c["alt"])] += 1

    # Runs that produced NO calls never appear above, and a caller that reports nothing on a
    # sample is exactly what a specificity column has to count. The run set is taken from the
    # widest caller, so a silent one is scored as silent rather than as absent.
    all_runs = collections.defaultdict(set)
    for c in calls:
        all_runs[(int(c["mapq_floor"]), int(c["min_reads"]), c["arm"])].add(c["run"])

    cols = ["caller", "mapq_floor", "min_reads", "arm", "certified_vaf", "runs_with_calls",
            "runs_expected", "calls_per_run_mean", "calls_per_run_max", "hotspot_detected",
            "hotspot_vaf_mean"]
    lines = ["\t".join(cols)]
    for key in sorted(burden, key=lambda k: (k[0], k[1], k[2], str(k[3]))):
        caller, mapq, mr, arm = key
        expected = len(all_runs[(mapq, mr, arm)])
        b = burden[key]
        h = hotspot[key]
        detected = sum(1 for x in h if x > 0)
        seen = [x for x in h if x > 0]
        lines.append("\t".join([
            caller, str(mapq), str(mr), arm,
            f"{ARMS.get(arm, float('nan')):g}",
            str(len(runs[key])), str(expected),
            f"{sum(b) / len(b):.2f}", str(max(b)),
            str(detected),
            f"{sum(seen) / len(seen):.6f}" if seen else "0",
        ]))
    out = "\n".join(lines)
    print(out)

    print("\n# substitution spectrum of the 0%-certified arm, all callers pooled", file=sys.stderr)
    total = sum(spectrum.values()) or 1
    for (ref, alt), n in spectrum.most_common(6):
        print(f"#   {ref}->{alt}\t{n}\t{100 * n / total:.1f}%", file=sys.stderr)

    if a.runtime and a.runtime.exists():
        secs = collections.defaultdict(list)
        for line in a.runtime.read_text().splitlines()[1:]:
            f = line.split("\t")
            if len(f) >= 6:
                secs[f[0]].append(float(f[4]))
        print("\n# wall clock per call, seconds", file=sys.stderr)
        for caller, v in sorted(secs.items()):
            v.sort()
            print(f"#   {caller}\tn={len(v)}\tmedian={v[len(v) // 2]:.1f}", file=sys.stderr)

    if a.out:
        a.out.write_text(out + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
