#!/usr/bin/env python3
"""Variant-vs-error ratios on the MIGEC spike-ins — the published test of whether UMI consensus
actually works.

2026-08-13.

Shugay et al. 2014 spiked three known IGH clones into a real B-cell library: EHEB, and two
deliberate variants at one and two substitutions from it. The question the paper poses is not "can
you remove errors" — anything can remove errors by discarding rare sequences — but

    can you remove PCR and sequencing error while KEEPING a real variant that is rarer than
    the error cloud around it?

The discriminating metric is the ratio of a real spike-in to the *worst error* at the same
substitution distance:

    Err1 = the most abundant junction exactly 1 substitution from EHEB, excluding V1
    Err2 = the most abundant junction exactly 2 substitutions from EHEB, excluding V2

Published values, and what any UMI pipeline has to reproduce:

    quantity     raw reads   standard processing   MIGEC (UMI consensus)
    EHEB/Err1        362           1041-1085              9007-24696
    V1/Err1         1.35             3.1-3.8                26.5-75.9
    V2/Err2         0.28             1.7-2.0                  4.6-6.2

V2 sits BELOW the worst 2-substitution error on raw reads (0.28x), so no abundance threshold can
separate them — that is the whole reason molecular barcodes exist. Moving V1/Err1 from ~1.4 to
26-76 is a change in the evidence, not in a threshold, and it is the acceptance gate for migec's
consensus assembly.

Usage:
    python scripts/spikein_ratio.py reads.fq.gz                    # baseline, raw reads
    python scripts/spikein_ratio.py consensus.fq.gz --label migec  # after assembly

The junction is located by its flanking anchors rather than by alignment, so this runs on raw
reads with no reference and no V/J calling.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

TRUTH = Path(__file__).resolve().parent.parent / "tests" / "data" / "migec_spikein_truth.tsv"


def load_truth(path: Path) -> dict[str, str]:
    with open(path) as fh:
        return {r["clone_id"]: r["junction_nt"] for r in csv.DictReader(fh, delimiter="\t")}


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return 10**9
    return sum(x != y for x, y in zip(a, b))


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def iter_reads(path: str):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.rstrip("\n")


def collect_junctions(path: str, parent: str, flank: int = 12) -> Counter:
    """Count junctions of the same length as `parent`, located by its conserved 3' anchor.

    ⚠ Anchor on the 3' end only. The obvious approach — require both the first and last `flank`
    bases of the parent — silently loses exactly what this script exists to measure: EHEB-V1
    differs from EHEB at position 4 and EHEB-V2 at positions 7-8, i.e. *inside* a 5' anchor. Both
    variants then report as zero and the metric looks perfect.

    Both orientations are searched: the Experiment 1 library carries this clone reverse
    complemented, and missing that also turns a real signal into a reported zero.
    """
    n = len(parent)
    right = parent[-flank:]
    right_rc = revcomp(parent[:flank])  # the 3' anchor of the reverse-complemented junction
    counts: Counter = Counter()
    for seq in iter_reads(path):
        for anchor, rc in ((right, False), (right_rc, True)):
            i = seq.find(anchor)
            if i < 0:
                continue
            end = i + len(anchor)
            start = end - n
            if start < 0:
                continue
            cand = seq[start:end]
            if len(cand) == n:
                counts[revcomp(cand) if rc else cand] += 1
            break
    return counts


def ratios(counts: Counter, truth: dict[str, str]) -> dict:
    parent = truth["EHEB"]
    v1, v2 = truth["EHEB-V1"], truth["EHEB-V2"]
    known = {parent, v1, v2}

    err1 = err2 = 0
    err1_seq = err2_seq = ""
    n_at = Counter()
    for seq, c in counts.items():
        d = hamming(seq, parent)
        if d:
            n_at[d] += 1
        if seq in known:
            continue
        if d == 1 and c > err1:
            err1, err1_seq = c, seq
        elif d == 2 and c > err2:
            err2, err2_seq = c, seq

    def ratio(a: int, b: int) -> float | None:
        return None if b == 0 else a / b

    return {
        "EHEB": counts.get(parent, 0),
        "EHEB-V1": counts.get(v1, 0),
        "EHEB-V2": counts.get(v2, 0),
        "Err1": err1,
        "Err1_seq": err1_seq,
        "Err2": err2,
        "Err2_seq": err2_seq,
        "distinct_junctions": len(counts),
        "n_at_1sub": n_at[1],
        "n_at_2sub": n_at[2],
        "n_at_3sub": n_at[3],
        "EHEB/Err1": ratio(counts.get(parent, 0), err1),
        "V1/Err1": ratio(counts.get(v1, 0), err1),
        "V2/Err2": ratio(counts.get(v2, 0), err2),
    }


# Published targets, for reporting only -- this script never decides pass/fail on its own.
PUBLISHED = {
    "raw": {"V1/Err1": (1.0, 4.0), "V2/Err2": (0.2, 2.0)},
    "migec": {"V1/Err1": (26.5, 75.9), "V2/Err2": (4.6, 6.2)},
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("reads", help="FASTQ(.gz): raw reads, or consensus reads after assembly")
    ap.add_argument("--label", default="raw", choices=sorted(PUBLISHED), help="which target band")
    ap.add_argument("--truth", default=str(TRUTH))
    ap.add_argument("--tsv", help="append a row to this TSV")
    args = ap.parse_args()

    truth = load_truth(Path(args.truth))
    counts = collect_junctions(args.reads, truth["EHEB"])
    if not counts:
        print(f"no junctions matching the EHEB flanks were found in {args.reads}", file=sys.stderr)
        print("check the orientation and that this library contains the spike-in", file=sys.stderr)
        return 1
    r = ratios(counts, truth)

    print(f"{args.reads}  ({r['distinct_junctions']} distinct junctions)")
    print(f"  EHEB          {r['EHEB']:>10,}")
    print(f"  EHEB-V1       {r['EHEB-V1']:>10,}   (1 substitution, real)")
    print(f"  EHEB-V2       {r['EHEB-V2']:>10,}   (2 substitutions, real)")
    print(f"  Err1          {r['Err1']:>10,}   worst error at 1 sub, of {r['n_at_1sub']}")
    print(f"  Err2          {r['Err2']:>10,}   worst error at 2 subs, of {r['n_at_2sub']}")
    print()
    band = PUBLISHED[args.label]
    for key in ("EHEB/Err1", "V1/Err1", "V2/Err2"):
        val = r[key]
        shown = "n/a (no competing error)" if val is None else f"{val:.2f}"
        target = band.get(key)
        note = ""
        if target and val is not None:
            note = "  within published band" if target[0] <= val <= target[1] else (
                f"  OUTSIDE published {target[0]}-{target[1]}"
            )
        print(f"  {key:<12} {shown:>12}{note}")

    if args.tsv:
        new = not Path(args.tsv).exists()
        with open(args.tsv, "a") as fh:
            if new:
                fh.write("file\tlabel\t" + "\t".join(k for k in r if not k.endswith("_seq")) + "\n")
            fh.write(
                f"{args.reads}\t{args.label}\t"
                + "\t".join(str(v) for k, v in r.items() if not k.endswith("_seq"))
                + "\n"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
