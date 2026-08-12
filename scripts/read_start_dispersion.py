#!/usr/bin/env python3
# 2026-08-13  X1: are the reads sharing a (cell barcode, UMI) co-terminal?
#
# Everything MIGEC does assumes they are. An amplicon MIG is one molecule amplified from fixed
# primers, so its reads start at the same base and a single ungapped consensus is the right model.
# A 3' GEX library is not that: the molecule is fragmented after capture, so reads sharing a UMI
# tile it at different offsets and may not overlap at all. If that is what the data shows, then a
# single consensus per UMI is meaningless for 10x and `assemble` must partition each group into
# overlap components first -- which is a design decision for M1, not an option to add later.
#
# Reads the BAM over https by byte range, so this costs a few hundred MB rather than the 4.8 GB
# the file weighs. Needs pysam (not a migec dependency -- this is a one-off experiment):
#
#     uv pip install pysam
#     python scripts/read_start_dispersion.py --bam <url-or-path> --region chr11:65497688-65508073

from __future__ import annotations

import argparse
import collections
import statistics
import sys


def collect(bam, regions, max_reads):
    """(cell, umi, contig, strand) -> [aligned blocks per read] over the given regions.

    Blocks, not (start, end): these are spliced RNA-seq alignments, so a read's genomic span
    includes its introns. Two reads either side of a junction can span the same megabase while
    sharing no aligned base at all, and `start - end` overlap would call them one contiguous
    fragment. ``get_blocks()`` is the M/=/X runs only.
    """
    import pysam

    groups: dict[tuple, list[list[tuple[int, int]]]] = collections.defaultdict(list)
    n_reads = n_tagged = 0
    with pysam.AlignmentFile(bam, "rb") as af:
        for region in regions:
            contig, _, span = region.partition(":")
            start, _, end = span.partition("-")
            for r in af.fetch(contig, int(start), int(end)):
                n_reads += 1
                # Secondary and supplementary alignments are the same fragment reported twice, and
                # a duplicate is the same fragment sequenced twice; either would fake co-terminal
                # reads that are really one read counted more than once.
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                    continue
                # CB/UB are the *corrected* barcodes. CR/UR are raw, and grouping on those would
                # measure barcode error as if it were fragmentation.
                if not r.has_tag("CB") or not r.has_tag("UB"):
                    continue
                n_tagged += 1
                # Strand matters: a molecule's reads all come off the same strand, and mixing them
                # would report a spurious spread between two genes on opposite strands.
                key = (r.get_tag("CB"), r.get_tag("UB"), r.reference_name, r.is_reverse)
                blocks = r.get_blocks()
                if blocks:
                    groups[key].append(blocks)
                if n_tagged >= max_reads:
                    return groups, n_reads, n_tagged
    return groups, n_reads, n_tagged


def shares_a_base(a, b):
    """Do two reads' aligned blocks intersect anywhere?"""
    for s1, e1 in a:
        for s2, e2 in b:
            if s1 < e2 and s2 < e1:
                return True
    return False


def overlap_components(reads):
    """How many groups of mutually-overlapping reads the molecule's reads fall into.

    This is the number the design turns on. One component means the reads cover one window and a
    single ungapped consensus is the right model -- MIGEC's assumption. More than one means the
    reads tile the molecule with gaps, and `assemble` has to emit one consensus per component or
    silently invent sequence across a gap it never observed.
    """
    parent = list(range(len(reads)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(reads)):
        for j in range(i + 1, len(reads)):
            if shares_a_base(reads[i], reads[j]):
                parent[find(i)] = find(j)
    return len({find(i) for i in range(len(reads))})


def summarise(groups, read_length):
    """One row per (CB,UMI) group with >1 read: how far apart its reads sit."""
    rows = []
    for blocks_per_read in groups.values():
        if len(blocks_per_read) < 2:
            continue
        starts = [b[0][0] for b in blocks_per_read]
        ends = [b[-1][1] for b in blocks_per_read]
        # Aligned bases only, so a shared intron does not count as shared sequence.
        covered = sorted({p for b in blocks_per_read for s, e in b for p in (s, e)})
        rows.append(
            {
                "reads": len(blocks_per_read),
                "distinct_starts": len(set(starts)),
                "start_range": max(starts) - min(starts),
                "end_range": max(ends) - min(ends),
                "coterminal": len(set(starts)) == 1,
                "components": overlap_components(blocks_per_read),
                # A footprint wider than one read means the group tiles the molecule rather than
                # covering one window of it -- the contig case, not the consensus case.
                "spans_gt_read": (covered[-1] - covered[0]) > read_length,
            }
        )
    return rows


def report(rows, n_reads, n_tagged, n_groups, read_length):
    multi = len(rows)
    if not multi:
        return "no (CB,UMI) group in this region carries more than one read"
    frac = lambda pred: sum(1 for r in rows if pred(r)) / multi  # noqa: E731
    reads_in_multi = sum(r["reads"] for r in rows)
    starts = [r["distinct_starts"] for r in rows]
    ranges = sorted(r["start_range"] for r in rows)
    comps = collections.Counter(r["components"] for r in rows)
    out = [
        f"reads examined            {n_reads:,}",
        f"  with CB and UB          {n_tagged:,}",
        f"(CB,UMI) groups           {n_groups:,}",
        f"  with >1 read            {multi:,} ({100 * multi / max(n_groups, 1):.1f}%), "
        f"holding {reads_in_multi:,} reads",
        "",
        f"reads per group           mean {statistics.mean(r['reads'] for r in rows):.2f}, "
        f"max {max(r['reads'] for r in rows)}",
        f"distinct starts per group mean {statistics.mean(starts):.2f}, max {max(starts)}",
        f"genomic start range (nt)  median {ranges[len(ranges) // 2]:,}, "
        f"p90 {ranges[int(0.9 * len(ranges))]:,}, max {ranges[-1]:,}",
        "",
        f"co-terminal groups                       {100 * frac(lambda r: r['coterminal']):5.1f}%",
        f"groups wider than one read ({read_length} nt)        "
        f"{100 * frac(lambda r: r['spans_gt_read']):5.1f}%",
        "",
        "overlap components per group (1 = a single ungapped consensus is valid):",
    ]
    for k in sorted(comps):
        out.append(f"  {k:>3}   {comps[k]:>7,}  {100 * comps[k] / multi:5.1f}%")
    out.append(f"  reads needing more than one consensus: "
               f"{100 * frac(lambda r: r['components'] > 1):.1f}% of groups")
    return "\n".join(out)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True, help="position-sorted BAM, path or https URL")
    p.add_argument(
        "--region",
        action="append",
        required=True,
        help="contig:start-end, repeatable. Pick expressed loci or nothing has >1 read per UMI.",
    )
    p.add_argument("--max-reads", type=int, default=2_000_000)
    p.add_argument("--read-length", type=int, default=91, help="10x v3 R2 is 91 nt")
    p.add_argument("--tsv", help="also write the per-group table here")
    a = p.parse_args(argv)

    groups, n_reads, n_tagged = collect(a.bam, a.region, a.max_reads)
    rows = summarise(groups, a.read_length)
    print(report(rows, n_reads, n_tagged, len(groups), a.read_length))

    if a.tsv:
        cols = ["reads", "distinct_starts", "start_range", "end_range", "coterminal",
                "components", "spans_gt_read"]
        with open(a.tsv, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(int(r[c]) if isinstance(r[c], bool) else r[c])
                                   for c in cols) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
