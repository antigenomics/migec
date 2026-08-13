#!/usr/bin/env python3
# 2026-08-13  Do the two models checkout reports actually predict what the reads show?
#
# checkout prints two numbers derived from models rather than from observation, and both feed
# decisions downstream:
#
#   1. the COLLISION rate, from the birthday problem on the barcode's collision entropy. It sets
#      `effective_space`, the saturation warning, and the corrected molecule count.
#   2. the UMI ERROR rate, estimated from the excess of barcode pairs at Hamming distance 1.
#
# Each can be checked against something independent. A collision is two molecules sharing a
# barcode, and if those molecules have different sequences the reads say so directly -- no model
# involved. A barcode error is either a sequencing miscall, which the reported Phred predicts, or
# a polymerase miscall, which the enzyme's error rate and the cycle count predict.
#
#     python scripts/collision_check.py --checkout out/CTRL.fq.gz --json out/checkout.json

from __future__ import annotations

import argparse
import collections
import gzip
import json
import math
import random
import sys


def read_tagged(path, window):
    """(umi, umi_quality, payload) per read, from checkout's RX/QX tags."""
    umi = qual = None
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                f = line.split()
                umi = next((x[5:] for x in f if x.startswith("RX:Z:")), None)
                qual = next((x[5:] for x in f if x.startswith("QX:Z:")), None)
            elif i % 4 == 1:
                s = line.rstrip("\n")
                if umi and len(s) >= window:
                    yield umi, qual, s[:window]


def differs(a, b, frac):
    """Do two sequences differ at more than `frac` of positions?"""
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] != b[i]) > frac * n


def split_by_template(reads, frac, min_share):
    """Greedy clustering of a MIG's reads into distinct templates.

    Greedy, not exhaustive: a collision is two sequences differing at many positions, which no
    reasonable clustering disagrees about. Only clusters holding `min_share` of the reads count,
    so a single erroneous read is not a second template.
    """
    reps = []
    for r in reads:
        for c in reps:
            if not differs(c[0], r, frac):
                c[1] += 1
                break
        else:
            reps.append([r, 1])
    return [c for c in reps if c[1] >= min_share * len(reads)]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--checkout", required=True, help="a per-sample FASTQ written by checkout")
    p.add_argument("--json", required=True, help="checkout.json from the same run")
    p.add_argument("--window", type=int, default=180)
    p.add_argument("--min-reads", type=int, default=4, help="MIGs below this cannot show a split")
    p.add_argument("--divergence", type=float, default=0.05)
    p.add_argument("--min-share", type=float, default=0.20)
    p.add_argument("--pcr-error", type=float, default=1e-5, help="per base per cycle")
    p.add_argument("--pcr-cycles", type=int, default=25)
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args(argv)

    sample = json.load(open(a.json))["samples"][0]
    L = sample["umi_length"]
    space = sample["effective_space"]
    l_eff = sample["effective_length"]

    migs = collections.defaultdict(list)
    quals = collections.Counter()
    for umi, qual, seq in read_tagged(a.checkout, a.window):
        migs[umi].append(seq)
        if qual:
            for ch in qual:
                quals[ord(ch) - 33] += 1
    n_obs = len(migs)

    # ---------------------------------------------------------------- collisions vs the birthday
    # Molecules land in the barcode space independently, so occupancy per barcode is Poisson. The
    # observed count of *distinct barcodes* is the occupied count, which is what pins lambda.
    #     occupied = S (1 - e^-lambda),  M = S lambda
    lam = -math.log1p(-n_obs / space) if n_obs < space else float("inf")
    m_hat = space * lam
    p_multi = (1 - math.exp(-lam) - lam * math.exp(-lam)) / (1 - math.exp(-lam))

    print("BARCODE SPACE")
    print(f"  barcode length            {L} nt  ({4**L:,} sequences)")
    print(f"  effective length          {l_eff:.3f} nt  ({space:,.0f} sequences)")
    print("    -- the collision entropy, not Shannon: sum_a p_a^2 per position, not -sum p log p")
    print(f"  distinct barcodes seen    {n_obs:,}")
    print(f"  occupancy                 {100 * n_obs / space:.1f}%")
    print(f"  Poisson lambda            {lam:.4f} molecules per barcode")
    print(f"  molecules implied         {m_hat:,.0f}   ({m_hat - n_obs:,.0f} hidden by collision)")
    print(f"\n  birthday prediction: {100 * p_multi:.2f}% of occupied barcodes hold >1 molecule")

    # Independent check: a collision of two *different* templates is visible in the reads.
    rng = random.Random(a.seed)
    big = {u: v for u, v in migs.items() if len(v) >= a.min_reads}
    split = sum(1 for v in big.values()
                if len(split_by_template(v, a.divergence, a.min_share)) > 1)
    # ...but only a collision between molecules that *differ* is visible, and in a clonal-ish
    # population most pairs do not. Estimate that detection probability from the data itself.
    reps = [v[0] for v in big.values()]
    trials = min(20_000, len(reps) * (len(reps) - 1) // 2)
    hits = sum(differs(rng.choice(reps), rng.choice(reps), a.divergence) for _ in range(trials))
    p_detect = hits / trials if trials else 0.0

    print("\n  observed, from the sequences (no model):")
    print(f"    MIGs with >= {a.min_reads} reads      {len(big):,}")
    print(f"    holding >1 template     {split:,}  ({100 * split / max(len(big), 1):.2f}%)")
    print(f"    P(two molecules differ) {p_detect:.4f}  -- measured on random pairs")
    if p_detect > 0:
        corrected = split / len(big) / p_detect
        print(f"    -> collision rate       {100 * corrected:.2f}%  after correcting for the "
              f"{100 * (1 - p_detect):.1f}% of collisions between identical templates")
        ratio = corrected / p_multi if p_multi else 0
        print(f"\n  observed / birthday       {ratio:.2f}x")

    # ------------------------------------------------------------------- barcode error vs models
    print("\nBARCODE ERROR RATE")
    tot_q = sum(quals.values())
    if tot_q:
        e_phred = sum(n * 10 ** (-q / 10) for q, n in quals.items()) / tot_q
        print(f"  from reported Phred       {e_phred:.3e} per base")
        print(f"    mean Q {sum(q * n for q, n in quals.items()) / tot_q:.1f}, "
              f"{len(quals)} distinct values over {tot_q:,} barcode bases")
        print("    -- the mean of 10^-Q/10, not 10^-(mean Q)/10; the low-Q tail dominates")
    e_pcr = a.pcr_error * a.pcr_cycles
    print(f"  from PCR                  {e_pcr:.3e} per base "
          f"({a.pcr_error:.0e} x {a.pcr_cycles} cycles)")
    if tot_q:
        print(f"  the two together          {e_phred + e_pcr:.3e} per base")
    print(f"  migec's estimate          {sample['umi_error_rate']:.3e} per base")
    print("    -- from the excess of barcode pairs at Hamming distance 1")
    if tot_q:
        pred = e_phred + e_pcr
        print(f"\n  estimated / predicted     {sample['umi_error_rate'] / pred:.2f}x")
        print(f"  expected barcodes with >=1 error: "
              f"{100 * (1 - (1 - pred) ** L):.2f}% at {L} nt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
