#!/usr/bin/env python3
# 2026-08-13  X2: fit the RT/PCR error floor on a clonal control.
#
# A consensus over c reads suppresses *sequencing* error like 1/c. It suppresses nothing that was
# already in the molecule when the first PCR cycle started: an RT miscall, or an early polymerase
# error, is in every read of the MIG and every consensus reproduces it faithfully. So the residual
# error of a consensus behaves like
#
#     e_out(c) = p_floor + a/c
#
# and the *intercept* is the floor -- the error rate no amount of over-sequencing removes. It is
# the single most consequential constant in the project, because it caps every quality migec emits
# above ~Q40, and it is currently a guess spanning 10x in both directions (1e-4 / 1e-5 / 1e-6).
# Measure it instead.
#
# The control is 8E5: a cell line carrying one integrated HIV-1 provirus, so every molecule in the
# library is the same sequence and any disagreement between a MIG's consensus and the library's
# modal sequence is error we introduced.
#
#     python scripts/quality_floor.py --reads SRR1763769_1.fastq.gz --out x2/

from __future__ import annotations

import argparse
import collections
import gzip
import math
import pathlib
import sys


def open_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fastq(path, limit=None):
    with open_maybe_gz(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                name = line[1:].split()[0]
            elif i % 4 == 1:
                seq = line.rstrip("\n")
            elif i % 4 == 3:
                yield name, seq, line.rstrip("\n")
                if limit and (i // 4) + 1 >= limit:
                    return


def cycle_profile(path, n_cycles, sample):
    """Per-cycle base composition, entropy and consensus base.

    This is what `migec suggest` will do: a UMI block reads as ~2 bits of entropy per cycle, a
    primer as ~0, and the boundary between them is the pattern. Reading it off the data beats
    trusting a protocol description, which is written for the bench and not for the file.
    """
    counts = [collections.Counter() for _ in range(n_cycles)]
    n = 0
    for _, seq, _ in read_fastq(path, sample):
        for j, b in enumerate(seq[:n_cycles]):
            counts[j][b] += 1
        n += 1
    prof = []
    for j, c in enumerate(counts):
        tot = sum(c.values()) or 1
        h = -sum((v / tot) * math.log2(v / tot) for v in c.values() if v)
        top, top_n = c.most_common(1)[0]
        prof.append({"cycle": j, "entropy": h, "consensus": top, "frac": top_n / tot})
    return prof, n


def infer_pattern(prof, umi_min_entropy=1.5, primer_min_frac=0.9, primer_len=20):
    """A leading run of high-entropy cycles, then the conserved primer that anchors it."""
    umi = 0
    while umi < len(prof) and prof[umi]["entropy"] >= umi_min_entropy:
        umi += 1
    if umi == 0:
        raise SystemExit(
            "no high-entropy block at the read start -- the UMI is not here. Check the other mate."
        )
    primer = []
    for p in prof[umi : umi + primer_len]:
        if p["frac"] < primer_min_frac:
            break
        primer.append(p["consensus"])
    if len(primer) < 8:
        raise SystemExit(
            f"only {len(primer)} conserved cycles after the {umi} nt UMI -- nothing to anchor on"
        )
    # Lowercase: a scored position at half weight, which is what the adapter region is.
    return "N" * umi + "".join(primer).lower(), umi, len(primer)


def consensus(seqs):
    """Majority base per column, over the shortest common length. No indel handling anywhere.

    A tie is ``N``, not an arbitrary pick. Ties are the whole story at even MIG sizes -- two reads
    disagreeing have no majority, and `Counter.most_common` would resolve it by insertion order,
    turning a coin flip into a confident wrong base. Scoring skips ``N``.
    """
    n = min(len(s) for s in seqs)
    out = []
    for j in range(n):
        c = collections.Counter(s[j] for s in seqs).most_common(2)
        out.append("N" if len(c) > 1 and c[0][1] == c[1][1] else c[0][0])
    return "".join(out)


def minor_allele_spectrum(consensuses, length):
    """Per position, the fraction of molecules not carrying the modal base.

    One vote per molecule, never per read: a single over-amplified molecule would otherwise carry
    its own RT error into the reference and hide it. ``N`` (an unresolved tie) does not vote.
    """
    out = []
    for j in range(length):
        c = collections.Counter(s[j] for s in consensuses if s[j] != "N")
        tot = sum(c.values())
        out.append(0.0 if not tot else 1.0 - c.most_common(1)[0][1] / tot)
    return out


def poisson_ci(k, n):
    """Exact-ish Poisson 95% interval for a rate k/n, as a rate. Wilson would do; this is simpler
    and correct in the regime that matters, which is k small or zero."""
    if n == 0:
        return 0.0, 0.0
    # Garwood: chi2 quantiles, approximated by Wilson-Hilferty so nothing extra is imported.
    def chi2(p, df):
        if df == 0:
            return 0.0
        z = 1.959963984540054 * (1 if p > 0.5 else -1)
        return df * (1 - 2 / (9 * df) + z * math.sqrt(2 / (9 * df))) ** 3
    lo = chi2(0.025, 2 * k) / 2 if k else 0.0
    hi = chi2(0.975, 2 * (k + 1)) / 2
    return lo / n, hi / n


def plateau_floor(bins, min_size):
    """Pool the residual over MIGs of at least `min_size` reads.

    Not a least-squares fit of ``p_floor + a/c``. That model is wrong for a majority-vote
    consensus: the sequencing residual is the chance that *most* reads carry the same wrong base,
    which falls roughly geometrically in c rather than as 1/c. Regressing on 1/c lets the smallest,
    noisiest bin set the intercept, and on simulated data with a known floor it returned a
    *negative* probability. The floor is simply where the curve flattens, so measure it there and
    give it an interval.
    """
    mm = sum(v[0] for k, v in bins.items() if k >= min_size)
    bases = sum(v[1] for k, v in bins.items() if k >= min_size)
    rate = mm / bases if bases else 0.0
    return rate, mm, bases, poisson_ci(mm, bases)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reads", required=True, help="clonal-control FASTQ carrying the UMI")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--window", type=int, default=120, help="payload bases to score per read")
    p.add_argument("--sample", type=int, default=200_000, help="reads used for the cycle profile")
    p.add_argument("--min-mig", type=int, default=2)
    p.add_argument("--max-mig", type=int, default=200, help="above this a MIG is a UMI collision")
    p.add_argument(
        "--max-minor",
        type=float,
        default=0.01,
        help="a position whose minor allele reaches this fraction of molecules is real variation, "
        "not error, and is excluded. 0.01 suits a viral quasispecies; a truly clonal control can "
        "use a smaller value.",
    )
    p.add_argument(
        "--max-divergence",
        type=float,
        default=0.05,
        help="a MIG disagreeing with the reference at more than this fraction of positions is a "
        "different template, not an erroneous copy, and is excluded. At any plausible floor the "
        "expected disagreement is well under 1%%, so this cut is unambiguous.",
    )
    a = p.parse_args(argv)

    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    prof, n_profiled = cycle_profile(a.reads, 80, a.sample)
    pattern, umi_len, primer_len = infer_pattern(prof)
    print(f"cycle profile over {n_profiled:,} reads")
    print(f"  UMI      {umi_len} nt at cycles 0-{umi_len - 1} "
          f"(mean entropy {sum(x['entropy'] for x in prof[:umi_len]) / umi_len:.2f} bits)")
    print(f"  primer   {primer_len} nt, {pattern[umi_len:].upper()}")
    print(f"  pattern  {pattern}")
    with open(out / "cycle_profile.tsv", "w") as fh:
        fh.write("cycle\tentropy_bits\tconsensus\tfraction\n")
        for x in prof:
            fh.write(f"{x['cycle']}\t{x['entropy']:.4f}\t{x['consensus']}\t{x['frac']:.4f}\n")

    # Extract with our own checkout rather than a bespoke parser: it is the code that has to be
    # right, and running it on real data is the point.
    from migec.checkout import run

    (out / "barcodes.txt").write_text(f"CTRL\t{pattern}\n")
    summary = run(a.reads, out / "barcodes.txt", out / "checkout")
    print(f"\ncheckout  {summary['assigned']:,} / {summary['total']:,} reads assigned "
          f"({100 * summary['assigned'] / summary['total']:.1f}%)")

    migs = collections.defaultdict(list)
    with gzip.open(out / "checkout" / "CTRL.fq.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                # The header is "@name<SP>RX:Z:...<TAB>QX:Z:...<TAB>BC:Z:...", so the name and the
                # first tag share a field if you split on tabs alone.
                umi = next((f[5:] for f in line.split() if f.startswith("RX:Z:")), None)
            elif i % 4 == 1:
                s = line.rstrip("\n")
                if umi and len(s) >= a.window:
                    migs[umi].append(s[: a.window])
    print(f"          {len(migs):,} distinct UMIs, "
          f"{sum(len(v) for v in migs.values()):,} reads at >= {a.window} nt")

    usable = {u: v for u, v in migs.items() if a.min_mig <= len(v) <= a.max_mig}
    cons = {u: consensus(v) for u, v in usable.items()}
    print(f"          {len(usable):,} MIGs with {a.min_mig}-{a.max_mig} reads")
    if len(usable) < 100:
        raise SystemExit("too few MIGs to fit anything")

    # Truth: the modal base per position across MIG *consensuses*, one vote per molecule. Voting
    # per read would let a single over-amplified molecule carry its RT error into the truth.
    truth = consensus(list(cons.values()))

    # ...but "disagrees with the modal base" is only error if the sample is clonal. An HIV plasma
    # population is a quasispecies: a position carrying a real 20% variant would contribute 0.2 to
    # the "error" rate and swamp a floor of 1e-4. So measure only where the molecules agree.
    #
    # A real variant below the threshold is still counted as error, which biases the floor *up* --
    # the safe direction for a quality cap, and it is stated rather than hidden.
    # The threshold has to sit well above 1/(number of molecules), or a position where a single
    # molecule carries an error looks "polymorphic" and is excluded -- which drops exactly the
    # positions the floor lives at and biases it *down*, the unsafe direction.
    if a.max_minor * len(cons) < 10:
        raise SystemExit(
            f"--max-minor {a.max_minor} over {len(cons)} molecules excludes any position with "
            f"{a.max_minor * len(cons):.1f} erroneous molecules, which is the signal itself. "
            f"Use at least {10 / len(cons):.3g}, or more molecules."
        )
    spectrum = minor_allele_spectrum(list(cons.values()), len(truth))
    monomorphic = [j for j, minor in enumerate(spectrum) if minor < a.max_minor]
    polymorphic = [j for j, minor in enumerate(spectrum) if minor >= a.max_minor]
    print(f"          {len(monomorphic)} of {len(truth)} positions monomorphic at "
          f"< {a.max_minor:.0%} minor allele; {len(polymorphic)} excluded as real variation")
    if polymorphic:
        worst = sorted(polymorphic, key=lambda j: -spectrum[j])[:8]
        print("          excluded positions: "
              + ", ".join(f"{j}({spectrum[j]:.1%})" for j in worst)
              + (" ..." if len(polymorphic) > 8 else ""))
    with open(out / "position_spectrum.tsv", "w") as fh:
        fh.write("position\ttruth_base\tminor_fraction\tmonomorphic\n")
        for j, minor in enumerate(spectrum):
            fh.write(f"{j}\t{truth[j]}\t{minor:.6f}\t{int(minor < a.max_minor)}\n")
    if len(monomorphic) < 20:
        raise SystemExit("fewer than 20 monomorphic positions -- this sample is not clonal enough")

    # A molecule that disagrees with the reference at a large *fraction* of positions is not an
    # erroneous copy of it, it is a different template -- another region, an off-target product,
    # or an indel-shifted read, and we model no indels anywhere. Over 180 bases at a floor of
    # 1e-4 the expected count is 0.02, so 5% (9 mismatches) is not a rate this process produces.
    #
    # This has to be cut, not averaged in: on the HIV control below, 0.8% of MIGs sit past 20%
    # divergence and contribute 90% of every mismatch in the dataset. Left in, they set the
    # "floor" two orders of magnitude too high and it is not an error rate at all.
    divergence = {
        u: sum(1 for j in monomorphic if c[j] != "N" and c[j] != truth[j]) / len(monomorphic)
        for u, c in cons.items()
    }
    print(f"\n{'divergence from the reference':<34}{'MIGs':>8}{'share':>8}{'of mismatches':>15}")
    tot_mm = sum(divergence.values()) or 1
    edges = [(0.0, "exact"), (0.01, "<= 1%"), (0.05, "1-5%"), (0.20, "5-20%"), (1.01, "> 20%")]
    prev = -1e-9
    for hi, label in edges:
        sel = [d for d in divergence.values() if prev < d <= hi] if hi else []
        if hi == 0.0:
            sel = [d for d in divergence.values() if d == 0.0]
        print(f"  {label:<32}{len(sel):>8}{100 * len(sel) / len(cons):>7.1f}%"
              f"{100 * sum(sel) / tot_mm:>14.1f}%")
        prev = hi
    kept = {u for u, d in divergence.items() if d <= a.max_divergence}
    print(f"\n  keeping {len(kept):,} of {len(cons):,} MIGs at <= {a.max_divergence:.0%} divergence "
          f"({len(cons) - len(kept):,} excluded as different templates)")
    cons = {u: c for u, c in cons.items() if u in kept}
    if len(cons) < 100:
        raise SystemExit("too few MIGs left after excluding divergent templates")

    # e_out(c): bases of a MIG's consensus that disagree with truth, binned by MIG size. Positions
    # where either side is an unresolved tie are not scored -- a coin flip is not an error rate.
    bins = collections.defaultdict(lambda: [0, 0])  # size -> [mismatches, bases]
    n_migs = collections.Counter()
    for u, c in cons.items():
        n = len(usable[u])
        n_migs[n] += 1
        for j in monomorphic:
            if c[j] == "N" or truth[j] == "N":
                continue
            bins[n][1] += 1
            bins[n][0] += c[j] != truth[j]

    with open(out / "quality_floor.tsv", "w") as fh:
        fh.write("mig_size\tmigs\tbases\tmismatches\te_out\te_lo\te_hi\n")
        for size in sorted(bins):
            mm, bases = bins[size]
            e = mm / bases if bases else 0.0
            lo, hi = poisson_ci(mm, bases)
            fh.write(f"{size}\t{n_migs[size]}\t{bases}\t{mm}\t{e:.6e}\t{lo:.6e}\t{hi:.6e}\n")

    print(f"\n{'MIG size':>9}{'MIGs':>8}{'bases':>11}{'mismatch':>10}{'e_out':>12}"
          f"{'95% CI':>25}")
    for size in sorted(bins):
        mm, bases = bins[size]
        if bases < 10_000:
            continue
        lo, hi = poisson_ci(mm, bases)
        print(f"{size:>9}{n_migs[size]:>8}{bases:>11,}{mm:>10}{mm / bases:>12.3e}"
              f"{f'[{lo:.2e}, {hi:.2e}]':>25}")

    # The floor is the plateau, so it needs a plateau to sit on. Below 3 reads there is no
    # majority to speak of and the residual is still dominated by sequencing error.
    print()
    for min_size in (3, 5, 10):
        rate, mm, bases, (lo, hi) = plateau_floor(bins, min_size)
        q = -10 * math.log10(hi) if hi > 0 else float("inf")
        print(f"  MIGs >= {min_size:>2} reads:  p_floor {rate:.3e}  "
              f"[{lo:.2e}, {hi:.2e}]  over {mm} mismatches in {bases:,} bases"
              f"   -> Q cap {q:.1f}")

    # How much of the answer is the divergence cut? If the floor moves with the threshold, the cut
    # is doing the measuring rather than removing contaminants, and the number is not a floor.
    print(f"\n  sensitivity to --max-divergence (currently {a.max_divergence:.0%}):")
    for thr in (0.02, 0.05, 0.10, 0.20):
        sub = collections.defaultdict(lambda: [0, 0])
        for u, c in cons.items():
            if divergence[u] > thr:
                continue
            n = len(usable[u])
            for j in monomorphic:
                if c[j] != "N" and truth[j] != "N":
                    sub[n][1] += 1
                    sub[n][0] += c[j] != truth[j]
        r, m, b, (lo_, hi_) = plateau_floor(sub, 5)
        print(f"    <= {thr:>4.0%}  p_floor {r:.3e}  [{lo_:.2e}, {hi_:.2e}]  ({m} in {b:,})")

    rate, mm, bases, (lo, hi) = plateau_floor(bins, 5)
    (out / "fit.tsv").write_text(
        f"p_floor\t{rate:.6e}\np_floor_lo\t{lo:.6e}\np_floor_hi\t{hi:.6e}\n"
        f"mismatches\t{mm}\nbases\t{bases}\n"
        f"q_cap\t{-10 * math.log10(hi) if hi > 0 else 0:.2f}\n"
    )
    print("\nThe Q cap is computed from the *upper* bound: claiming a quality the data cannot")
    print("support is the failure that matters, and at these counts the point estimate is loose.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
