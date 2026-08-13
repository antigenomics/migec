#!/usr/bin/env python3
# 2026-08-13  Everything a UMI library will not tell you unless asked.
#
# Five questions, in the order you actually want them answered:
#
#   1. WHAT IS THIS?  Align the assembled consensus to a reference and say which gene and which
#      coordinates it covers. A pipeline that reports 124,501 molecules without saying what they
#      are is reporting a number, not a result.
#   2. IS THE UMI ANY GOOD?  The per-position PWM of the barcode itself -- not the per-cycle PWM
#      `suggest` uses to FIND it, but the composition of the extracted barcodes, which is what the
#      collision arithmetic runs on.
#   3. WHERE ON THE FLOWCELL?  Illumina headers carry lane, tile and cluster coordinates. Reads of
#      one UMI sitting on top of each other are one cluster read twice -- an optical duplicate --
#      which inflates the MIG size and fakes over-sequencing.
#   4. IS THE INDEX HOPPING?  A Casava 1.8 header carries the sample index. The i7 x i5 table is
#      the only way hopping is estimable, because a hopped read carries a *valid* index and an
#      "ambiguous" bucket therefore measures nothing.
#   5. WHAT VARIES?  The per-position minor-allele spectrum of the consensuses against their own
#      modal sequence, which separates a quasispecies from an error process.
#
#     python scripts/diagnose.py --consensus asm/CTRL.consensus.fq.gz \
#         --checkout co/CTRL.fq.gz --reference HXB2.fasta

from __future__ import annotations

import argparse
import collections
import gzip
import math
import pathlib
import re
import sys

BASES = "ACGT"

# Casava 1.8+:  instrument:run:flowcell:lane:tile:x:y  [space] read:filtered:control:index
CASAVA18 = re.compile(
    r"^(?P<instrument>[^:]+):(?P<run>\d+):(?P<flowcell>[^:]+):(?P<lane>\d+):"
    r"(?P<tile>\d+):(?P<x>\d+):(?P<y>\d+)$"
)
# Casava <1.8:  instrument:lane:tile:x:y#index/read
CASAVA14 = re.compile(
    r"^(?P<instrument>[^:]+):(?P<lane>\d+):(?P<tile>\d+):(?P<x>\d+):(?P<y>-?\d+)"
    r"(?:#(?P<index>[^/]*))?(?:/(?P<read>\d))?$"
)


def parse_header(name, comment):
    """Lane/tile/x/y/index from an Illumina read name, or None.

    SRA normalises headers to `@SRR1763769.1 1/2` and the coordinates are gone for good -- so this
    returns None rather than guessing, and the caller reports that the diagnostic is unavailable
    instead of quietly reporting zero optical duplicates.
    """
    m = CASAVA18.match(name)
    if m:
        d = m.groupdict()
        fields = comment.split(":") if comment else []
        d["index"] = fields[3] if len(fields) > 3 else None
        return d
    m = CASAVA14.match(name)
    if m:
        return m.groupdict()
    return None


def read_fastq(path, want_tags=()):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        name = comment = seq = None
        for i, line in enumerate(fh):
            if i % 4 == 0:
                head = line[1:].rstrip("\n")
                name, _, comment = head.partition(" ")
            elif i % 4 == 1:
                seq = line.rstrip("\n")
            elif i % 4 == 3:
                tags = {}
                for f in comment.split("\t"):
                    for sep in (":Z:", ":i:"):
                        if sep in f:
                            k, v = f.split(sep, 1)
                            if not want_tags or k in want_tags:
                                tags[k] = v
                yield name, comment, seq, line.rstrip("\n"), tags


# ------------------------------------------------------------------ 1. what is this


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def load_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                chunks.append(line.strip())
    return name, "".join(chunks).upper()


def best_ungapped(query, reference):
    """Best ungapped placement of `query` in `reference`, either strand.

    Ungapped on purpose: migec models no indels anywhere, and a gapped aligner here would paper
    over exactly the frameshifts that matter. Returns (mismatches, start, strand, aligned length).
    """
    best = None
    for strand, probe in (("+", query), ("-", revcomp(query))):
        n = len(probe)
        for i in range(len(reference) - n + 1):
            d = 0
            for a, b in zip(probe, reference[i : i + n]):
                if a != b and a != "N" and b != "N":
                    d += 1
                    if best and d >= best[0]:
                        break
            else:
                if best is None or d < best[0]:
                    best = (d, i, strand, n)
    return best


HXB2_GENES = [
    ("5' LTR", 1, 634), ("gag", 790, 2292), ("pol", 2085, 5096), ("protease", 2253, 2549),
    ("RT", 2550, 3869), ("RNase H", 3870, 4229), ("integrase", 4230, 5096), ("vif", 5041, 5619),
    ("vpr", 5559, 5850), ("tat", 5831, 6045), ("vpu", 6062, 6310), ("env", 6225, 8795),
    ("nef", 8797, 9417), ("3' LTR", 9086, 9719),
]


def identify(consensus, reference, ref_name):
    hit = best_ungapped(consensus, reference)
    if hit is None:
        return ["no placement found"]
    d, start, strand, n = hit
    a, b = start + 1, start + n
    out = [
        f"reference   {ref_name[:70]}",
        f"placement   {a:,}-{b:,} ({n} nt) on the {strand} strand, "
        f"{d} mismatches ({d / n:.1%} divergence)",
    ]
    overlapping = [g for g, s, e in HXB2_GENES if not (b < s or a > e)]
    if overlapping:
        out.append(f"covers      {', '.join(overlapping)}")
    return out


# ------------------------------------------------------------------ 2. the UMI's own PWM


def umi_pwm(umis):
    if not umis:
        return []
    length = len(next(iter(umis)))
    rows = []
    for j in range(length):
        c = collections.Counter(u[j] for u in umis if u[j] in BASES)
        n = sum(c.values()) or 1
        freq = {b: c[b] / n for b in BASES}
        h = -sum(f * math.log2(f) for f in freq.values() if f)
        m = sum(f * f for f in freq.values())
        rows.append({"position": j, **freq, "entropy": h, "collision": m,
                     "effective_bases": 1 / m if m else 0.0})
    return rows


# ------------------------------------------------------- 3/4. flowcell and index hopping


def flowcell_report(records, radius):
    """Optical duplicates, per-tile yield, and the i7 x i5 table when the header carries one."""
    coords = collections.defaultdict(list)   # (umi) -> [(lane, tile, x, y)]
    per_tile = collections.Counter()
    indices = collections.Counter()
    parsed = total = 0
    for name, comment, _seq, _qual, tags in records:
        total += 1
        h = parse_header(name, comment)
        if not h:
            continue
        parsed += 1
        per_tile[(h["lane"], h["tile"])] += 1
        if h.get("index"):
            indices[h["index"]] += 1
        umi = tags.get("RX")
        if umi:
            coords[umi].append((h["lane"], h["tile"], int(h["x"]), int(h["y"])))
    out = {"reads": total, "with_coordinates": parsed, "tiles": len(per_tile),
           "indices": indices, "optical_duplicates": 0, "groups_checked": 0}
    for umi, pts in coords.items():
        if len(pts) < 2:
            continue
        out["groups_checked"] += 1
        by_tile = collections.defaultdict(list)
        for lane, tile, x, y in pts:
            by_tile[(lane, tile)].append((x, y))
        for pts_here in by_tile.values():
            for i in range(len(pts_here)):
                for j in range(i + 1, len(pts_here)):
                    dx = pts_here[i][0] - pts_here[j][0]
                    dy = pts_here[i][1] - pts_here[j][1]
                    if dx * dx + dy * dy <= radius * radius:
                        out["optical_duplicates"] += 1
                        break
    return out


def index_table(indices):
    """i7 x i5 contingency. A hopped read carries a *valid* index, so the signal is the
    combinations that were never declared -- an 'ambiguous' bucket measures nothing at all."""
    pairs = collections.Counter()
    for idx, n in indices.items():
        i7, _, i5 = idx.partition("+")
        pairs[(i7, i5 or ".")] += n
    i7s = sorted({a for a, _ in pairs})
    i5s = sorted({b for _, b in pairs})
    return pairs, i7s, i5s


# ------------------------------------------------------------------ 5. what varies


def minor_spectrum(seqs, width):
    out = []
    for j in range(width):
        c = collections.Counter(s[j] for s in seqs if j < len(s) and s[j] in BASES)
        n = sum(c.values())
        if not n:
            out.append((0.0, "N", 0))
            continue
        top, top_n = c.most_common(1)[0]
        out.append((1 - top_n / n, top, n))
    return out


# ------------------------------------------------------------------------ driver


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--consensus", required=True, help="assemble's consensus FASTQ")
    p.add_argument("--checkout", help="checkout's tagged FASTQ, for the flowcell diagnostics")
    p.add_argument("--reference", help="FASTA to place the consensus against")
    p.add_argument("--optical-radius", type=int, default=100,
                   help="pixels within which two clusters of one UMI are one cluster read twice")
    p.add_argument("--max-consensus", type=int, default=20000,
                   help="Consensus records to use. Sampled uniformly across the file, NEVER the "
                        "first N: assemble writes in barcode order, so the first N all begin with "
                        "the same letters and the barcode PWM comes out fixed at its leading "
                        "positions. Measured: the first 4,000 of this HIV library report an "
                        "effective length of 6.45 nt against the true 8.97.")
    p.add_argument("--out")
    a = p.parse_args(argv)

    # Reservoir sample, because assemble writes in barcode order and the first N records of a
    # sorted file share their leading bases. Taking them would report a barcode that is fixed at
    # its first positions -- 6.45 effective nt against the true 8.97 on the library below -- which
    # is the same mistake as subsampling reads instead of whole barcodes, wearing a third hat.
    import random

    rng = random.Random(0)
    reservoir, seen = [], 0
    for _name, _c, seq, _q, tags in read_fastq(a.consensus):
        seen += 1
        item = (seq, tags.get("RX"), int(tags["cD"]) if "cD" in tags else None)
        if len(reservoir) < a.max_consensus:
            reservoir.append(item)
        else:
            j = rng.randrange(seen)
            if j < a.max_consensus:
                reservoir[j] = item
    if not reservoir:
        raise SystemExit("no consensus records read")
    cons = [r[0] for r in reservoir]
    umis = [r[1] for r in reservoir if r[1]]
    depths = [r[2] for r in reservoir if r[2] is not None]
    print(f"(sampled {len(reservoir):,} of {seen:,} consensus records, uniformly)")

    # The width most molecules reach, not the shortest one: a single short consensus would
    # otherwise truncate the whole diagnostic to its own length, and on a library that is 56%
    # singletons the shortest is very short indeed.
    lengths = sorted(len(s) for s in cons)
    width = lengths[len(lengths) // 4]     # the length three quarters of molecules reach
    spectrum = minor_spectrum(cons, width)
    modal = "".join(t for _, t, _ in spectrum)

    print("=" * 78)
    print("WHAT IS THIS SAMPLE?")
    print("=" * 78)
    print(f"consensuses {len(cons):,} sampled; lengths {lengths[0]}-{lengths[-1]} nt, "
          f"median {lengths[len(lengths) // 2]}, scoring the {width} nt that three quarters reach")
    print(f"modal       {modal[:70]}{'...' if len(modal) > 70 else ''}")
    if a.reference:
        ref_name, ref = load_fasta(a.reference)
        for line in identify(modal, ref, ref_name):
            print(line)
    else:
        print("(no --reference given, so the sample is not identified)")

    print()
    print("=" * 78)
    print("THE UMI'S OWN PWM")
    print("=" * 78)
    rows = umi_pwm(umis)
    if rows:
        print(f"{'pos':>4}{'A':>8}{'C':>8}{'G':>8}{'T':>8}{'bits':>8}{'eff bases':>11}")
        for r in rows:
            print(f"{r['position']:>4}{r['A']:>8.3f}{r['C']:>8.3f}{r['G']:>8.3f}{r['T']:>8.3f}"
                  f"{r['entropy']:>8.3f}{r['effective_bases']:>11.3f}")
        total_bits = sum(r["entropy"] for r in rows)
        eff_len = -sum(math.log(r["collision"], 4) for r in rows if r["collision"] > 0)
        print(f"\n{'':>4}{'':>32}{total_bits:>8.3f}  bits total")
        print(f"effective length {eff_len:.3f} nt of {len(rows)} -- "
              f"{4 ** eff_len:,.0f} usable barcodes of {4 ** len(rows):,}")
    else:
        print("(no RX tags on the consensus records)")

    print()
    print("=" * 78)
    print("THE FLOWCELL")
    print("=" * 78)
    if not a.checkout:
        print("(no --checkout given)")
    else:
        fc = flowcell_report(read_fastq(a.checkout, want_tags=("RX",)), a.optical_radius)
        print(f"reads       {fc['reads']:,}")
        if not fc["with_coordinates"]:
            print("coordinates NONE. This file's headers carry no lane/tile/x/y -- SRA normalises")
            print("            them to `@SRR....N N/2` and they are gone for good. Optical")
            print("            duplicates and index hopping are BOTH unmeasurable here; that is a")
            print("            property of the file, not a clean result.")
        else:
            print(f"coordinates {fc['with_coordinates']:,} reads "
                  f"({fc['with_coordinates'] / fc['reads']:.1%}) over {fc['tiles']} tiles")
            print(f"optical     {fc['optical_duplicates']:,} of {fc['groups_checked']:,} "
                  f"multi-read UMIs hold two clusters within {a.optical_radius} px")
            print("            -- one cluster read twice inflates the MIG and fakes")
            print("            over-sequencing, so it is not a molecule seen twice")
            if fc["indices"]:
                pairs, i7s, i5s = index_table(fc["indices"])
                print(f"\ni7 x i5     {len(i7s)} x {len(i5s)} = {len(i7s) * len(i5s)} combinations, "
                      f"{len(pairs)} seen")
                print(f"{'':>12}{' '.join(f'{b[:8]:>10}' for b in i5s)}")
                for i7 in i7s:
                    cells = " ".join(f"{pairs.get((i7, i5), 0):>10,}" for i5 in i5s)
                    print(f"{i7[:10]:>12}{cells}")
                declared = sum(n for (x, y), n in pairs.items() if n > 0.01 * max(pairs.values()))
                hopped = sum(fc["indices"].values()) - declared
                print(f"\nunused combinations carry {hopped:,} reads -- the only estimable")
                print("signal, because a hopped read carries a VALID index and an 'ambiguous'")
                print("bucket therefore measures nothing")
            else:
                print("index       the header carries none, so hopping is not estimable")

    print()
    print("=" * 78)
    print("WHAT VARIES")
    print("=" * 78)
    polymorphic = [(j, m, t) for j, (m, t, _n) in enumerate(spectrum) if m >= 0.01]
    print(f"{len(polymorphic)} of {width} positions carry a minor allele at >= 1% of molecules")
    # Where they sit decides what they are. Real variation is scattered; read-end quality decay
    # piles up in the last cycles, and reporting a "quasispecies" that is really the 3' tail of the
    # read is the mistake this line exists to prevent.
    if polymorphic:
        tail = [j for j, _m, _t in polymorphic if j >= 0.85 * width]
        head_rate = sum(m for j, m, _t in polymorphic if j < 0.85 * width) / max(
            len([1 for j, _m, _t in polymorphic if j < 0.85 * width]), 1)
        tail_rate = sum(m for j, m, _t in polymorphic if j >= 0.85 * width) / max(len(tail), 1)
        print(f"  first 85% of the read: mean minor allele {head_rate:.2%}")
        print(f"  last  15% of the read: mean minor allele {tail_rate:.2%}"
              + (f"  <- {tail_rate / head_rate:.1f}x higher, which is read-end quality decay "
                 f"rather than biology" if head_rate and tail_rate > 2 * head_rate else ""))
    if polymorphic:
        print(f"\n{'pos':>5}{'modal':>7}{'minor':>9}")
        for j, m, t in sorted(polymorphic, key=lambda r: -r[1])[:15]:
            print(f"{j:>5}{t:>7}{m:>8.1%}")
    if depths:
        depths.sort()
        print(f"\nMIG depth   median {depths[len(depths) // 2]}, "
              f"p90 {depths[int(0.9 * len(depths))]}, max {depths[-1]}")

    if a.out:
        out = pathlib.Path(a.out)
        out.mkdir(parents=True, exist_ok=True)
        with open(out / "umi_pwm.tsv", "w") as fh:
            fh.write("position\tA\tC\tG\tT\tentropy_bits\tcollision\teffective_bases\n")
            for r in rows:
                fh.write(f"{r['position']}\t{r['A']:.6f}\t{r['C']:.6f}\t{r['G']:.6f}\t"
                         f"{r['T']:.6f}\t{r['entropy']:.6f}\t{r['collision']:.6f}\t"
                         f"{r['effective_bases']:.6f}\n")
        with open(out / "variation.tsv", "w") as fh:
            fh.write("position\tmodal\tminor_fraction\tmolecules\n")
            for j, (m, t, n) in enumerate(spectrum):
                fh.write(f"{j}\t{t}\t{m:.6f}\t{n}\n")
        print(f"\nwrote {out}/umi_pwm.tsv and {out}/variation.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
