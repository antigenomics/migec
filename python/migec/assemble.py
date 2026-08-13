"""The assemble stage: reads sharing a UMI are reads of one molecule, so collapse them.

All of the work is in C++ (`_core.assemble`); this module calls it once and turns the returned
summary into tables and a report. The per-molecule table is written by the C++ side because it is
one line per molecule and there are hundreds of millions of them.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from migec import _core
from migec.checkout import _bytes, _dur, _pct


# The pre-amplification error floor, by what made it. An error in the RT step or the first PCR
# cycle is in every read of the molecule, so no consensus removes it and no emitted quality may
# exceed -10 log10 of it. Which of these applies is a property of the PROTOCOL, not of the data,
# which is why it is named rather than fitted.
#
#   rt      1e-4   a reverse transcription step. 10x states it for V(D)J: "The estimated error rate
#                  for the V(D)J RT reaction is 1e-4 per base", and X2 measured 1.54e-4 on an
#                  HIV-1 Primer ID library independently. Illumina quote 7.37e-5 for TruSight
#                  Oncology 500 v2, the same decade.
#   medium  1e-5   no RT, an ordinary polymerase. McInerney et al. 2014 measured Taq at
#                  4.3e-5 +/- 1.8 per bp per template duplication; only the first cycle is
#                  common-mode to the molecule, so the decade below it.
#   high    1e-6   no RT, a proofreading polymerase. Same paper: Pfu 2.8e-6, Phusion 2.6e-6,
#                  Pwo 2.4e-6 per bp per duplication.
#
# Anything else: give the number. Never guess between them -- a floor set an order of magnitude
# too low is a Q60 printed on a base that is wrong one time in ten thousand.
#
# Note: the first cycle is not an ordinary cycle. Shagin et al. 2017 (Sci Rep 7:2718) measured
# nine polymerases at 0.3-6.6e-5 per base per cycle AND found that errors from the initial linear
# amplification are 5 +/- 1 times more frequent than the per-cycle rate of the PCR that follows.
# That first cycle is exactly what this floor is: it is copied into every read of the molecule,
# and the classes above are per-protocol brackets around it, not a fitted value. `--rt-error auto`
# (M1, open) is what will fit it per dataset. See docs/quality_floor.rst and SOURCES.md.
RT_FLOORS = {"rt": 1e-4, "medium": 1e-5, "high": 1e-6}


def parse_rt_error(value: str) -> float:
    """`rt` / `medium` / `high`, or a number. Raises ValueError naming all four."""
    key = value.strip().lower()
    if key in RT_FLOORS:
        return RT_FLOORS[key]
    try:
        floor = float(key)
    except ValueError:
        raise ValueError(
            f"--rt-error {value!r} is neither a fidelity class nor a number. Give one of "
            + ", ".join(f"{n} ({v:.0e}, caps at Q{-10 * math.log10(v):.0f})"
                        for n, v in RT_FLOORS.items())
            + ", or the rate itself"
        ) from None
    if not 0.0 < floor < 1.0:
        raise ValueError(f"--rt-error {value!r} is not a per-base error rate in (0, 1)")
    return floor


def _quantiles_by_depth(grid: list[dict]) -> list[tuple]:
    """Five-number summary of emitted quality, per power-of-two depth bin.

    The grid is the exact joint count of (depth bin, rounded Phred), so these are real order
    statistics over every molecule -- not a summary of a thinned sample. Quality is discrete and
    capped at the RT floor, which is why a box is the honest shape here and a scatter is not: at
    high depth every molecule sits on the same one or two integers, and a cloud of dots draws that
    as a line whether it holds ten molecules or ten million.
    """
    bins: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for row in grid:
        bins.setdefault((row["min_reads"], row["max_reads"]), []).append(
            (row["quality"], row["molecules"])
        )

    def at(counts: list[tuple[int, int]], total: int, fraction: float) -> int:
        # Nearest-rank. The value returned is always a quality that molecules actually had.
        target, seen = max(1, math.ceil(fraction * total)), 0
        for quality, n in counts:
            seen += n
            if seen >= target:
                return quality
        return counts[-1][0]

    out = []
    for (lo, hi), counts in sorted(bins.items()):
        counts.sort()
        total = sum(n for _, n in counts)
        mean = sum(q * n for q, n in counts) / total
        out.append(
            (
                lo,
                hi,
                total,
                counts[0][0],
                at(counts, total, 0.25),
                at(counts, total, 0.50),
                at(counts, total, 0.75),
                counts[-1][0],
                f"{mean:.3f}",
            )
        )
    return out


def run(
    reads: str | Path,
    out_dir: str | Path,
    sample_id: str = "",
    # Both are measured constants and both live in consensus.hpp, next to what measured them.
    rt_floor: float = _core.RT_FLOOR,
    linkage_threshold: float = _core.LINKAGE_THRESHOLD,
    contig: bool = False,
    fast: bool = False,
    min_reads: int = 1,
    gzip_level: int = _core.GZIP_LEVEL,
    threads: int = 0,
    limit_reads: int = 0,
    limit_umis: int = 0,
) -> dict:
    """Assemble consensuses from a checkout-tagged FASTQ into `out_dir`."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary = _core.assemble(
        str(reads), str(out), sample_id, rt_floor, linkage_threshold, contig, fast, min_reads,
        gzip_level, 0, threads, limit_reads, limit_umis,
    )
    summary["input"] = str(reads)
    summary["rt_floor"] = rt_floor
    summary["linkage_threshold"] = linkage_threshold

    with open(out / "assemble.coverage.tsv", "w") as fh:
        fh.write("sample_id\tmin_reads\tmax_reads\tgroups\n")
        for b in summary["coverage"]:
            fh.write(
                f"{summary['sample_id']}\t{b['min_reads']}\t{b['max_reads']}\t{b['groups']}\n"
            )

    with open(out / "assemble.quality_by_depth.tsv", "w") as fh:
        fh.write(
            "sample_id\tmin_reads\tmax_reads\tmolecules\tq_min\tq_p25\tq_median\tq_p75\tq_max\t"
            "q_mean\n"
        )
        for row in _quantiles_by_depth(summary["quality_grid"]):
            fh.write(f"{summary['sample_id']}\t" + "\t".join(str(v) for v in row) + "\n")
    (out / "assemble.json").write_text(json.dumps(summary, indent=2, default=str))
    return summary


def format_report(summary: dict) -> str:
    s = summary
    lines = [
        f"reads       {s['reads']:,}",
        f"  grouped   {s['reads'] - s['reads_without_umi'] - s['reads_dropped']:,} into "
        f"{s['groups']:,} UMIs",
    ]
    if s["reads_without_umi"]:
        lines.append(
            f"  no RX tag  {s['reads_without_umi']:,} "
            f"({_pct(s['reads_without_umi'], s['reads'])}) -- was this file written by checkout?"
        )
    if s["reads_dropped"]:
        lines.append(f"  dropped   {s['reads_dropped']:,} in groups below --min-reads")
    if s["groups_capped"]:
        lines.append(
            f"  capped    {s['reads_over_cap']:,} reads in {s['groups_capped']:,} groups over "
            f"{s['max_reads_per_group']:,} reads did not enter the consensus -- they are still "
            f"counted as reads of their molecule"
        )
    lines += [
        "",
        f"molecules   {s['molecules']:,}",
        f"  split     {s['groups_split']:,} groups held more than one "
        f"({_pct(s['groups_split'], max(s['groups'], 1))} of groups)",
        f"  expected  {s['expected_molecules_per_group']:.2f} molecules per group from the "
        f"birthday problem at {s['barcode_space']['occupancy']:.1%} occupancy, "
        f"{s['molecules'] / max(s['groups'], 1):.2f} recovered",
        "",
        f"{_dur(s['wall_seconds'])} = {_dur(s['partition_seconds'])} partitioning into "
        f"{s['buckets']} bucket(s), serial + "
        f"{_dur(s['wall_seconds'] - s['partition_seconds'])} consensus on {s['threads']} threads",
        f"peak RSS {_bytes(s['peak_rss_bytes'])} -- one bucket is resident at a time",
        "",
        f"mean emitted quality  Q{s['mean_quality']:.1f}  "
        f"(capped at Q{s['quality_cap']:.0f} by the RT floor of {s['rt_floor']:.1e})",
        f"mean consensus error  {s['mean_consensus_error']:.2e} before the floor is added",
    ]
    if s["fast_mode"]:
        lines += [
            f"modal sequence carried {s['mean_support']:.1%} of the reads in the groups that "
            f"emitted one -- counting mode, so no per-base correction and no sub-clustering",
        ]
    lines += ["", f"{'MIG size':>12}{'groups':>12}{'share':>9}"]
    if s["contig_mode"]:
        lines.insert(
            len(lines) - 1,
            f"contigs     {s['contigs']:,} from {s['groups_fragmented']:,} groups whose reads did "
            f"not all reach each other "
            f"({_pct(s['groups_fragmented'], max(s['groups'], 1))} of groups)",
        )
    total = sum(b["groups"] for b in s["coverage"]) or 1
    for b in s["coverage"]:
        if not b["groups"]:
            continue
        label = f"{b['min_reads']}" if b["min_reads"] == b["max_reads"] else (
            f"{b['min_reads']}-{b['max_reads']}"
        )
        lines.append(f"{label:>12}{b['groups']:>12,}{100 * b['groups'] / total:>8.1f}%")

    warnings = []
    if s.get("limited"):
        warnings.append(
            "the intake was limited (--limit-read / --limit-umi), so every number here describes "
            "the FIRST reads of the file rather than the library. Nothing measured under a limit "
            "-- error rate, occupancy, molecule count -- transfers to the whole run"
        )
    singletons = next((b["groups"] for b in s["coverage"] if b["min_reads"] == 1), 0)
    if singletons > 0.5 * total:
        warnings.append(
            f"{_pct(singletons, total)} of molecules were seen once. A consensus over one read is "
            f"that read -- the UMI is buying counting here, not error correction"
        )
    # Splitting is the decision with the worst failure mode: a split molecule is counted twice and
    # nothing downstream can tell. The threshold is a measured false-positive point, so a rate far
    # above it means the threshold is being applied to data it was not measured on.
    if s["groups"] and s["groups_split"] / s["groups"] > 0.05:
        warnings.append(
            f"{_pct(s['groups_split'], s['groups'])} of groups were split into more than one "
            f"molecule. X3 measured a 1.60% rate on a real library at this threshold "
            f"({s['linkage_threshold']:.2f}); much more than that usually means UMI collisions, "
            f"not subclones -- check `p_multi` in checkout's barcode space table"
        )
    # A short UMI cannot tag every input molecule distinctly, and that is a design choice rather
    # than a defect -- but it decides whether a consensus per barcode means anything, and it makes
    # contig assembly unsafe, because two fragments of two different molecules on one barcode
    # share no sequence and look exactly like two fragments of one.
    if s["expected_molecules_per_group"] > 1.1:
        bound = "at least " if s["barcode_space"]["saturated"] else ""
        warnings.append(
            f"a group holds {bound}{s['expected_molecules_per_group']:.2f} molecules on average at "
            f"{s['barcode_space']['occupancy']:.0%} occupancy of a "
            f"{s['barcode_space']['effective_length']:.1f} nt effective barcode. The consensus of "
            f"a group is then a mixture of templates, not a molecule -- resolve them by sequence "
            f"(the split above) or use a longer UMI"
        )
    if s["contig_mode"] and s["barcode_space"]["p_multi"] > 0.05:
        warnings.append(
            f"contig assembly at {s['barcode_space']['p_multi']:.0%} of groups holding more than "
            f"one molecule: two fragments of two different molecules on one barcode have no "
            f"sequence in common, which is indistinguishable from two fragments of one. The "
            f"contigs here are not trustworthy"
        )
    if warnings:
        lines.append("")
        lines += [f"warning: {w}" for w in warnings]
    return "\n".join(lines)
