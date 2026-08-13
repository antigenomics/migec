#!/usr/bin/env python3
# 2026-08-13  X3: replace three derivations with three permutations.
#
# Three numbers in this pipeline come out of a derivation that assumes something the data has
# never been asked about. Each one has a permutation that measures the same quantity while
# assuming nothing, and the difference between the two is the size of the assumption.
#
#   A. `P_coll = prod_j sum_a p_j(a)^2` assumes the barcode's POSITIONS ARE INDEPENDENT. A
#      synthesiser mixes each coupling step separately, so they need not be. Null: compare the
#      collision probability of the observed joint distribution over k adjacent positions against
#      the product of its own marginals, for k = 1, 2, 3, ...  The excess per added position
#      extrapolates to full length.
#
#   B. The barcode error rate is read off the EXCESS of pairs at Hamming distance 1 over a chance
#      background, and both the background and the excess are derived. Two nulls: shuffle the
#      columns (marginals kept, error children destroyed) to measure the chance background, and
#      shuffle the READ COUNTS over the fixed distance-1 graph (graph kept, parent/child count
#      relationship destroyed) to measure how many of those pairs are really parent and child.
#
#   C. Splitting a MIG into two consensuses is accepted at a threshold derived from a Poisson
#      argument. That argument treats reads as exchangeable, and they are not -- a bad read carries
#      a minor base at many positions at once, which looks exactly like a linked subclone. Null:
#      randomise the minor-allele matrix keeping BOTH margins (per-position error count and
#      per-read error count), and read the threshold off the false-positive curve.
#
#     python scripts/permutation_nulls.py --reads SRR1763769_1.fastq.gz --out x3/

from __future__ import annotations

import argparse
import collections
import gzip
import itertools
import json
import math
import pathlib
import random
import sys

BASES = "ACGT"


# --------------------------------------------------------------------------------- input


def load_migs(fastq, window):
    """UMI -> list of payload prefixes, from checkout's RX tag."""
    migs = collections.defaultdict(list)
    umi = None
    with gzip.open(fastq, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                # "@name<SP>RX:Z:...<TAB>QX:Z:..." -- splitting on tabs alone merges name and RX.
                umi = next((f[5:] for f in line.split() if f.startswith("RX:Z:")), None)
            elif i % 4 == 1:
                s = line.rstrip("\n")
                if umi and len(s) >= window:
                    migs[umi].append(s[:window])
    return migs


# ------------------------------------------------------------- A. the independence null


def collision(counts):
    """The probability two independent draws coincide, as an UNBIASED estimate.

    `sum_a p_hat(a)^2` is biased upward by `(1 - sum p^2) / n`, and that bias grows as the
    distribution spreads -- so it grows with the k-mer width, and reads as dependence that
    accumulates with k. On the HIV library it produced a fake 1.007x over 9 nt out of nothing.
    `sum_a n_a(n_a - 1) / (N(N - 1))` is the U-statistic and has no such bias.

    Callers pass barcodes with N already folded to A, exactly as `pack_barcode` stores them: the
    packed key is what every stage groups on, so it is the alphabet the collision rate is about.
    Counting N as a fifth letter instead lets `m_j` fall below 1/4, which is impossible for a
    four-letter distribution -- the bound is the check, and it shows up as an effective length
    longer than the barcode.
    """
    values = list(counts.values())
    n = sum(values)
    return sum(v * (v - 1) for v in values) / (n * (n - 1)) if n > 1 else 0.0


def jensen_shannon(counts, q, total):
    """JSD in bits between an observed k-mer distribution and the product of its own marginals.

    The independence null *is* a distribution -- q(u) = prod_j p_j(u_j), a product measure -- so
    the honest test compares the whole of P against the whole of Q, not one functional of each.
    KL(P || Q) over the full barcode is exactly the multi-information; JSD is its bounded
    symmetric form and does not blow up on a cell the null gives zero weight.
    """
    total = float(total)
    out = 0.0
    for word, mass in q.items():
        p = counts.get(word, 0) / total
        m = 0.5 * (p + mass)
        if p > 0:
            out += 0.5 * p * math.log2(p / m)
        if mass > 0:
            out += 0.5 * mass * math.log2(mass / m)
    return out


def product_measure(marginals):
    """q(u) = prod_j p_j(u_j) over every word of the window: the independence null itself."""
    q = {"": 1.0}
    for freq in marginals:
        q = {w + b: m * freq[b] for w, m in q.items() for b in BASES if freq[b] > 0}
    return q


def jsd_null(barcodes, L, rng, reps, widths):
    """Observed JSD to the product measure, against the JSD a sample of the same size gets by
    chance. The floor matters: on 4^k cells and n barcodes the empirical distribution is sparse,
    so JSD is bounded away from zero even under perfect independence. A column shuffle draws from
    the product measure at exactly the observed n, which is that floor, measured rather than
    derived.
    """
    n = len(barcodes)
    freqs = []
    for j in range(L):
        c = collections.Counter(b[j] for b in barcodes)
        freqs.append({b: c[b] / n for b in BASES})
    shuffles = [column_shuffle(barcodes, rng) for _ in range(reps)]
    rows = []
    for k in widths:
        obs, null = [], []
        for j in range(0, L - k + 1):
            q = product_measure(freqs[j : j + k])
            obs.append(jensen_shannon(
                collections.Counter(b[j : j + k] for b in barcodes), q, n))
            for sh in shuffles:
                null.append(jensen_shannon(
                    collections.Counter(b[j : j + k] for b in sh), q, len(sh)))
        mu = sum(null) / len(null)
        sd = math.sqrt(sum((x - mu) ** 2 for x in null) / max(len(null) - 1, 1))
        mean_obs = sum(obs) / len(obs)
        rows.append({"k": k, "windows": len(obs), "jsd_observed": mean_obs,
                     "jsd_null_mean": mu, "jsd_null_sd": sd,
                     "excess": mean_obs - mu, "z": (mean_obs - mu) / sd if sd else 0.0})
    return rows


def independence_null(barcodes, L, max_k=None):
    """Collision of the observed joint over k adjacent positions, against the product of marginals.

    Both sides are computed on the same barcode list, so the comparison is internal and needs no
    model. It runs on *distinct* barcodes, which is the conservative direction: a barcode drawn
    twice appears once, so the high-frequency combinations that dependence creates are undercounted
    and the excess reported here is a lower bound.
    """
    n = len(barcodes)
    marg = [collision(collections.Counter(b[j] for b in barcodes)) for j in range(L)]
    # A k-mer distribution needs many more barcodes than the 4^k cells it spreads over, or its
    # collision is dominated by sampling noise rather than by dependence.
    if max_k is None:
        max_k = 1
        while 4 ** (max_k + 1) <= n / 100 and max_k < L:
            max_k += 1
    rows = []
    for k in range(1, max_k + 1):
        ratios = []
        for j in range(L - k + 1):
            obs = collision(collections.Counter(b[j : j + k] for b in barcodes))
            ind = math.prod(marg[j : j + k])
            if ind > 0 and obs > 0:
                ratios.append(math.log10(obs / ind))
        rows.append({"k": k, "windows": len(ratios),
                     "log10_excess": sum(ratios) / len(ratios) if ratios else 0.0})
    # Excess accumulates per *added* position, so the slope over k is the per-position cost of
    # assuming independence and it extrapolates to the full barcode.
    slope = sum(r["log10_excess"] / (r["k"] - 1) for r in rows[1:]) / max(len(rows) - 1, 1)
    return {
        "marginal_collision": marg,
        "rows": rows,
        "log10_excess_per_position": slope,
        "predicted_excess_full_length": 10 ** (slope * (L - 1)),
        "p_coll_independent": math.prod(marg),
        "p_coll_dependent": math.prod(marg) * 10 ** (slope * (L - 1)),
    }


# ------------------------------------------------------- B. the distance-1 graph nulls


def d1_pairs(seqs):
    """Unordered pairs of sequences at Hamming distance 1, by 3L neighbour enumeration.

    Enumeration, not all-pairs: 3L lookups per barcode against n^2/2 comparisons. Index order
    breaks the double count -- each pair is emitted once, by its lower-indexed member.
    """
    index = {s: i for i, s in enumerate(seqs)}
    out = []
    for i, s in enumerate(seqs):
        for j in range(len(s)):
            for b in BASES:
                if b == s[j]:
                    continue
                k = index.get(s[:j] + b + s[j + 1 :])
                if k is not None and k > i:
                    out.append((i, k))
    return out


def column_shuffle(barcodes, rng):
    """New barcodes with every position's marginal preserved and all dependence destroyed."""
    cols = [[b[j] for b in barcodes] for j in range(len(barcodes[0]))]
    for c in cols:
        rng.shuffle(c)
    return list({"".join(t) for t in zip(*cols)})


def solve_epsilon(counts, excess, L, lo=1e-7, hi=0.2):
    """Invert  excess = 3L * sum_i [ t_i + t_i^2 ],  t_i = 1 - exp(-c_i * eps / 3).

    eps/3, not eps: a miscall lands on one specific one of the three other bases, so a given
    neighbour of the parent is reached at a third of the per-base rate. The squared term is the
    sibling pair -- two children of one parent are themselves at distance 1.
    """
    def expected(eps):
        tot = 0.0
        for c in counts:
            t = 1.0 - math.exp(-c * eps / 3.0)
            tot += t + t * t
        return 3 * L * tot

    if excess <= 0 or expected(hi) < excess:
        return 0.0 if excess <= 0 else hi
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if expected(mid) < excess:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def graph_nulls(barcodes, counts, L, rng, reps, ratios):
    """Chance background from a column shuffle; parent/child structure from a count shuffle."""
    obs = d1_pairs(barcodes)

    # How many distance-1 pairs would this many barcodes of this composition have by chance?
    null_counts = []
    for _ in range(reps):
        shuffled = column_shuffle(barcodes, rng)
        # The shuffle can collide, so scale to the observed barcode count: pairs go as n^2.
        scale = (len(barcodes) / len(shuffled)) ** 2
        null_counts.append(len(d1_pairs(shuffled)) * scale)
    chance = sum(null_counts) / len(null_counts)

    # Of the pairs that exist, how many are parent and child? An error child is much smaller than
    # its parent; two unrelated neighbours are not. Permuting the counts over the *same* graph
    # keeps the composition, the graph and the count distribution, and destroys only that link.
    cvec = [counts[b] for b in barcodes]
    table = []
    for r in ratios:
        seen = sum(1 for i, j in obs if max(cvec[i], cvec[j]) >= r * min(cvec[i], cvec[j]))
        null = []
        for _ in range(reps):
            perm = cvec[:]
            rng.shuffle(perm)
            null.append(sum(1 for i, j in obs
                            if max(perm[i], perm[j]) >= r * min(perm[i], perm[j])))
        mu = sum(null) / len(null)
        sd = math.sqrt(sum((x - mu) ** 2 for x in null) / max(len(null) - 1, 1))
        table.append({"ratio": r, "observed": seen, "null_mean": mu, "null_sd": sd,
                      "excess": seen - mu, "z": (seen - mu) / sd if sd else float("inf")})

    return {
        "pairs_observed": len(obs),
        "pairs_by_chance": chance,
        "pairs_excess": len(obs) - chance,
        "epsilon_permutation": solve_epsilon(cvec, len(obs) - chance, L),
        "epsilon_all_pairs": solve_epsilon(cvec, len(obs), L),
        "count_ratio": table,
    }


# ------------------------------------------------------ C. the within-MIG linkage null


def log_choose(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(c11, n, a, b):
    """-log10 P(X >= c11) for X ~ Hypergeometric(n, a, b): a and b minors, c11 of them shared."""
    tot = -math.inf
    denom = log_choose(n, a)
    for x in range(c11, min(a, b) + 1):
        term = log_choose(b, x) + log_choose(n - b, a - x) - denom
        tot = term if tot == -math.inf else max(tot, term) + math.log1p(
            math.exp(-abs(tot - term))
        )
    return -tot / math.log(10) if tot > -math.inf else math.inf


def minor_matrix(reads, min_minor, max_cols):
    """Binary reads x positions matrix of "carries the minor base here"."""
    n = len(reads)
    width = min(len(r) for r in reads)
    cols = []
    for j in range(width):
        c = collections.Counter(r[j] for r in reads)
        major = c.most_common(1)[0][0]
        minor = n - c[major]
        if min_minor <= minor <= n // 2:
            cols.append((minor, j, major))
    cols.sort(reverse=True)
    cols = cols[:max_cols]
    return [[1 if r[j] != major else 0 for _, j, major in cols] for r in reads]


def linkage_score(mat):
    """Strongest co-segregation of minor alleles over any pair of columns, Bonferroni'd.

    Two-sided. At a 50/50 split which allele is the "major" one is a coin toss taken separately
    per column, so two columns of a genuine doublet can come out anti-correlated and a one-sided
    test scores the strongest evidence available as nothing. Minor-with-major is the same evidence
    in the opposite phase; both are tested and the count is halved for it.
    """
    n, k = len(mat), len(mat[0]) if mat else 0
    if k < 2:
        return 0.0
    colsum = [sum(row[c] for row in mat) for c in range(k)]
    best = 0.0
    correction = math.log10(k * (k - 1) // 2) + math.log10(2)
    for x, y in itertools.combinations(range(k), 2):
        c11 = sum(1 for row in mat if row[x] and row[y])
        a, b = colsum[x], colsum[y]
        same = hypergeom_sf(c11, n, a, b) if c11 >= 2 else 0.0
        opposite = hypergeom_sf(a - c11, n, a, n - b) if a - c11 >= 2 else 0.0
        best = max(best, max(same, opposite) - correction)
    return best


def curveball(mat, rng, swaps):
    """Randomise a binary matrix keeping BOTH margins (Strona et al. 2014, trade fashion).

    Both margins is the whole point. Permuting each column independently keeps the per-position
    error count but lets every read carry an average number of errors, and real reads do not --
    a low-quality read is minor at many positions at once and mimics a linked subclone exactly.
    """
    rows = [set(i for i, v in enumerate(r) if v) for r in mat]
    n = len(rows)
    if n < 2:
        return mat
    for _ in range(swaps):
        i, j = rng.randrange(n), rng.randrange(n)
        if i == j:
            continue
        shared = rows[i] & rows[j]
        a, b = rows[i] - shared, rows[j] - shared
        if not a or not b:
            continue
        pool = list(a | b)
        rng.shuffle(pool)
        rows[i] = shared | set(pool[: len(a)])
        rows[j] = shared | set(pool[len(a) :])
    k = len(mat[0])
    return [[1 if c in r else 0 for c in range(k)] for r in rows]


def linkage_null(migs, rng, min_reads, max_reads, min_minor, max_cols, reps, max_migs):
    keys = [u for u, v in migs.items() if min_reads <= len(v) <= max_reads]
    rng.shuffle(keys)
    keys = keys[:max_migs]
    obs, null = [], []
    for u in keys:
        mat = minor_matrix(migs[u], min_minor, max_cols)
        if not mat or len(mat[0]) < 2:
            continue
        obs.append(linkage_score(mat))
        swaps = 20 * len(mat)
        null += [linkage_score(curveball(mat, rng, swaps)) for _ in range(reps)]
    return keys, obs, null


def quantile(xs, q):
    if not xs:
        return 0.0
    s = sorted(xs)
    return s[min(len(s) - 1, max(0, int(round(q * (len(s) - 1)))))]


# ------------------------------------------------------------------------------ driver


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reads", required=True, help="FASTQ carrying the UMI")
    p.add_argument("--out", required=True)
    p.add_argument("--pattern", help="barcode pattern; inferred with `migec suggest` if omitted")
    p.add_argument("--cycles", type=int, default=40, help="cycles profiled by `migec suggest`")
    p.add_argument("--window", type=int, default=120, help="payload bases scored per read")
    p.add_argument("--reps", type=int, default=20, help="permutations per null")
    p.add_argument("--min-reads", type=int, default=10, help="MIG size for the linkage null")
    p.add_argument("--max-reads", type=int, default=200)
    p.add_argument("--min-minor", type=int, default=2, help="reads needed to call a position")
    p.add_argument("--max-cols", type=int, default=8)
    p.add_argument("--max-migs", type=int, default=2000)
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args(argv)

    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(a.seed)

    from migec.checkout import run
    from migec.suggest import run as suggest_run

    pattern = a.pattern
    if pattern is None:
        pattern = suggest_run(a.reads, out, cycles=a.cycles)["pattern"]
        print(f"suggest   {pattern}")
    (out / "barcodes.txt").write_text(f"CTRL\t{pattern}\n")
    summary = run(a.reads, out / "barcodes.txt", out / "checkout")
    sample = summary["samples"][0]
    L = sample["umi_length"]
    print(f"checkout  {summary['assigned']:,} / {summary['total']:,} reads assigned "
          f"({100 * summary['assigned'] / summary['total']:.1f}%), {L} nt UMI")

    migs = load_migs(out / "checkout" / "CTRL.fq.gz", a.window)
    # N -> A, as pack_barcode does. Every stage groups on the packed key, so that is the barcode
    # whose collision rate matters; the flag that a base was uncalled travels separately.
    fold = str.maketrans({c: "A" for c in "NnRYSWKMBDHV"})
    folded = collections.defaultdict(list)
    for u, v in migs.items():
        folded[u.translate(fold)] += v  # merge, never overwrite: N->A can land on a real barcode
    migs = folded
    barcodes = sorted(migs)
    counts = {u: len(v) for u, v in migs.items()}
    print(f"          {len(barcodes):,} distinct UMIs over "
          f"{sum(counts.values()):,} reads of >= {a.window} nt")

    # ------------------------------------------------------------------------------ A
    print("\nA. ARE THE BARCODE POSITIONS INDEPENDENT?")
    # Distinct barcodes, never barcodes weighted by reads: read counts are set by amplification,
    # so weighting would let one over-amplified molecule write the composition.
    #
    # ...but a hard read threshold throws away real molecules, which is the mistake this project
    # refuses everywhere else. Singletons are the compromise that costs least: an error child is
    # usually a singleton, so dropping them removes most of the barcodes that are not molecules
    # while keeping every molecule seen twice. Both are reported, and if they disagree the answer
    # is contamination by error children rather than dependence.
    subsets = [("all distinct", barcodes),
               ("reads >= 2", [b for b in barcodes if counts[b] >= 2])]
    ind = None
    for label, subset in subsets:
        if len(subset) < 5000:
            continue
        res = independence_null(subset, L)
        if ind is None:
            ind = res
        print(f"\n  {label}: {len(subset):,} barcodes")
        print(f"{'k':>6}{'windows':>9}{'observed/independent':>22}")
        for r in res["rows"]:
            print(f"{r['k']:>6}{r['windows']:>9}{10 ** r['log10_excess']:>21.5f}x")
        print(f"    per added position {10 ** res['log10_excess_per_position']:.5f}x, "
              f"over {L} positions {res['predicted_excess_full_length']:.4f}x")
        print(f"    effective length   {-math.log(res['p_coll_independent'], 4):.4f} nt "
              f"-> {-math.log(res['p_coll_dependent'], 4):.4f} nt under the null")

    # The collision rate is one functional of the distribution. The null IS a distribution -- the
    # product measure q(u) = prod_j p_j(u_j) -- so compare the whole of P against the whole of Q.
    print("\n  Jensen-Shannon divergence to the product measure (bits), against the floor a")
    print("  sample of the same size gets by chance:")
    jsd = jsd_null(barcodes, L, rng, max(a.reps // 5, 2), (2, 3, 4, 5))
    print(f"{'k':>6}{'observed':>12}{'chance':>12}{'sd':>10}{'excess':>12}{'z':>8}")
    for r in jsd:
        print(f"{r['k']:>6}{r['jsd_observed']:>12.3e}{r['jsd_null_mean']:>12.3e}"
              f"{r['jsd_null_sd']:>10.1e}{r['excess']:>12.3e}{r['z']:>8.1f}")
    print("  measured on distinct barcodes, so repeated draws are collapsed: a lower bound.")

    # ------------------------------------------------------------------------------ B
    print("\nB. HOW MANY DISTANCE-1 PAIRS ARE PARENT AND CHILD?")
    g = graph_nulls(barcodes, counts, L, rng, a.reps, (2, 5, 10, 20, 50))
    print(f"  pairs at distance 1   {g['pairs_observed']:,}")
    print(f"  expected by chance    {g['pairs_by_chance']:,.0f}   "
          f"(column shuffle, {a.reps} reps: marginals kept, error children destroyed)")
    print(f"  excess                {g['pairs_excess']:,.0f}")
    print(f"\n{'count ratio >=':>16}{'observed':>11}{'null mean':>12}{'null sd':>9}"
          f"{'excess':>11}{'z':>9}")
    for r in g["count_ratio"]:
        print(f"{r['ratio']:>16}{r['observed']:>11,}{r['null_mean']:>12,.0f}{r['null_sd']:>9.0f}"
              f"{r['excess']:>11,.0f}{r['z']:>9.1f}")
    print(f"\n  barcode error from the excess   {g['epsilon_permutation']:.3e} per base")
    print(f"  ...if every pair were a child   {g['epsilon_all_pairs']:.3e}  (the upper bound)")
    print(f"  checkout's estimate             {sample['error_budget']['estimated']:.3e}")
    print(f"  Phred + polymerase prediction   {sample['error_budget']['predicted']:.3e}")

    # ------------------------------------------------------------------------------ C
    print("\nC. WHEN IS A MIG REALLY TWO MOLECULES?")
    keys, obs, null = linkage_null(migs, rng, a.min_reads, a.max_reads, a.min_minor,
                                   a.max_cols, max(a.reps // 4, 1), a.max_migs)
    if not obs:
        print("  no MIG had two callable positions -- nothing to calibrate")
        thresholds = []
    else:
        print(f"  {len(obs):,} MIGs of {a.min_reads}-{a.max_reads} reads with >= 2 callable "
              f"positions, {len(null):,} randomisations")
        print(f"\n{'target FP':>11}{'threshold':>12}{'MIGs called':>13}{'called %':>10}")
        thresholds = []
        for fp in (0.05, 0.01, 0.001):
            t = quantile(null, 1 - fp)
            called = sum(1 for s in obs if s > t)
            thresholds.append({"target_fp": fp, "threshold": t, "called": called,
                               "fraction": called / len(obs)})
            print(f"{fp:>11.3f}{t:>12.2f}{called:>13,}{100 * called / len(obs):>9.2f}%")
        naive = sum(1 for s in obs if s > 2.0)  # -log10 p = 2 is the usual alpha = 0.01
        print(f"\n  a nominal p < 0.01 threshold (score > 2.00) would call {naive:,} "
              f"({100 * naive / len(obs):.2f}%)")
        print(f"  the permutation puts the 1% false-positive point at "
              f"{quantile(null, 0.99):.2f} instead")

    (out / "nulls.json").write_text(json.dumps(
        {"pattern": pattern, "umi_length": L, "distinct_umis": len(barcodes),
         "independence": ind, "jsd": jsd, "graph": g,
         "linkage": {"migs": len(obs), "randomisations": len(null),
                     "thresholds": thresholds,
                     "null_quantiles": {str(q): quantile(null, q)
                                        for q in (0.5, 0.9, 0.99, 0.999)}}},
        indent=2))
    with open(out / "nulls.linkage.tsv", "w") as fh:
        fh.write("source\tscore\n")
        for s in obs:
            fh.write(f"observed\t{s:.4f}\n")
        for s in null:
            fh.write(f"permuted\t{s:.4f}\n")
    print(f"\nwrote {out / 'nulls.json'} and {out / 'nulls.linkage.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
