#!/usr/bin/env python3
# 2026-08-14
#
# Call sets over the certified ctDNA arms, scored the same way whoever produced them.
#
# Three questions, one table, because they are all "how many calls, and which of them are real":
#
#   1. caller against caller on the IDENTICAL consensus BAM -- LoFreq against Mutect2. If both
#      report the same positions on the 0%-certified arm then the artifact is in the input rather
#      than in one caller's threshold, which is the thing a single caller cannot establish.
#   2. the remedies `docs/detection.rst` lists for the adapter read-through -- a MAPQ floor and
#      adapter trimming -- against `--min-reads`, which is the third.
#   3. PIPELINE against pipeline: UMIErrorCorrect replaces `assemble` rather than following it, so
#      it is scored from its own consensus and its rows say so.
#
# The 0% arm is the number that matters. Its variant frequency is zero by construction, so every
# call there is that pipeline's own false-positive rate on real chemistry rather than on a
# simulation.
#
# Never: score SUBSTITUTIONS against substitutions. migec emits no indels by design, and 56% of
# UMIErrorCorrect's PASS calls are deletions -- comparing raw call counts would score a class we
# do not produce. Both totals are in the table; the comparison is the substitution column.
#
# Usage:
#     python scripts/compare_callers.py \
#         --variants migec=callers_variants.tsv \
#         --variants migec-trimmed=trimming_variants.tsv \
#         --variants umierrorcorrect=uec_variants.tsv \
#         --runtime callers_runtime.tsv --out assets/ctdna_callers.tsv

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ctdna_persite import ARMS, read_tsv  # noqa: E402

# The certified hotspot: PIK3CA H1047R, the one locus whose frequency the vendor certifies.
HOTSPOT = ("3", "179234297")
BASES = set("ACGT")


def load(label: str, path: Path) -> list[dict]:
    """One source's calls, with the columns it does not carry filled in by what it means.

    A table with no `caller` column came from LoFreq; one with no `min_reads` came from a pipeline
    whose group-size cutoff is its own default. Both are named here rather than left blank, so a
    row can be read without knowing which file it came from.
    """
    out = []
    for r in read_tsv(path):
        if r.get("filtered", "PASS") not in ("PASS", "."):
            continue
        r["pipeline"] = label
        # A table with no `caller` column came from a pipeline that calls for itself, so the label
        # names it: `migec+lofreq` means that caller, a bare `umierrorcorrect` means its own.
        r.setdefault("caller", label.split("+")[-1])
        r.setdefault("mapq_floor", "0")
        r.setdefault("min_reads", "3")
        r["substitution"] = r["ref"] in BASES and r["alt"] in BASES
        out.append(r)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", action="append", required=True, metavar="LABEL=PATH",
                    help="a call table; repeat for each pipeline")
    ap.add_argument("--runtime", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--hotspot", default=f"{HOTSPOT[0]}:{HOTSPOT[1]}")
    a = ap.parse_args(argv)

    chrom, pos = a.hotspot.split(":")
    calls = []
    for spec in a.variants:
        label, _, path = spec.partition("=")
        if not path:
            raise SystemExit(f"--variants wants LABEL=PATH, got {spec!r}")
        calls += load(label, Path(path))
    if not calls:
        raise SystemExit("no PASS calls in any input")

    per_run: dict = collections.defaultdict(list)
    for c in calls:
        key = (c["pipeline"], c["caller"], int(c["mapq_floor"]), int(c["min_reads"]), c["arm"])
        per_run[(key, c["run"])].append(c)

    rows = collections.defaultdict(lambda: {"runs": set(), "all": [], "subs": [], "hs": []})
    spectrum: dict = collections.defaultdict(collections.Counter)
    for (key, run), cs in per_run.items():
        agg = rows[key]
        agg["runs"].add(run)
        agg["all"].append(len(cs))
        agg["subs"].append(sum(1 for c in cs if c["substitution"]))
        hit = [c for c in cs if c["chrom"] == chrom and c["pos"] == pos]
        agg["hs"].append(float(hit[0]["af"]) if hit else 0.0)
        if key[4] == "WT":
            for c in cs:
                if c["substitution"]:
                    spectrum[key[0]][(c["ref"], c["alt"])] += 1

    cols = ["pipeline", "caller", "mapq_floor", "min_reads", "arm", "certified_vaf", "runs",
            "calls_per_run", "substitutions_per_run", "hotspot_detected", "hotspot_vaf_mean"]
    lines = ["\t".join(cols)]
    for key in sorted(rows, key=lambda k: (k[0], k[1], k[2], k[3], str(k[4]))):
        pipeline, caller, mapq, mr, arm = key
        agg = rows[key]
        seen = [x for x in agg["hs"] if x > 0]
        lines.append("\t".join([
            pipeline, caller, str(mapq), str(mr), arm,
            f"{ARMS.get(arm, float('nan')):g}",
            str(len(agg["runs"])),
            f"{sum(agg['all']) / len(agg['all']):.2f}",
            f"{sum(agg['subs']) / len(agg['subs']):.2f}",
            str(len(seen)),
            f"{sum(seen) / len(seen):.6f}" if seen else "0",
        ]))
    out = "\n".join(lines)
    print(out)

    for pipeline, sp in sorted(spectrum.items()):
        total = sum(sp.values()) or 1
        to_g = sum(n for (_ref, alt), n in sp.items() if alt == "G")
        print(f"\n# {pipeline}: substitution spectrum of the 0%-certified arm, "
              f"{to_g / total:.0%} of it -> G", file=sys.stderr)
        for (ref, alt), n in sp.most_common(5):
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
