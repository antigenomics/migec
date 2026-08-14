"""The assemble stage: reads sharing a UMI are reads of one molecule, so collapse them.

All of the work is in C++ (`_core.assemble`); this module calls it once and turns the returned
summary into tables and a report. The per-molecule table is written by the C++ side because it is
one line per molecule and there are hundreds of millions of them.
"""

from __future__ import annotations

import collections
import gzip
import json
import math
import shutil
from pathlib import Path

from migec import _core
from migec.buckets import mig_buckets as _mig_buckets
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


def parse_rt_error(value: str) -> float | str:
    """`auto`, `rt` / `medium` / `high`, or a number. Raises ValueError naming all of them."""
    key = value.strip().lower()
    if key == "auto":
        return "auto"
    if key in RT_FLOORS:
        return RT_FLOORS[key]
    try:
        floor = float(key)
    except ValueError:
        raise ValueError(
            f"--rt-error {value!r} is neither a fidelity class nor a number. Give 'auto', one of "
            + ", ".join(f"{n} ({v:.0e}, caps at Q{-10 * math.log10(v):.0f})"
                        for n, v in RT_FLOORS.items())
            + ", or the rate itself"
        ) from None
    if not 0.0 < floor < 1.0:
        raise ValueError(f"--rt-error {value!r} is not a per-base error rate in (0, 1)")
    return floor


def poisson_ci(k: int, n: int) -> tuple[float, float]:
    """Exact-ish Poisson 95% interval for a rate k/n, as a rate.

    Garwood's chi-square interval with the quantiles taken from Wilson-Hilferty, so nothing extra
    is imported. Correct in the regime that matters here, which is k small or zero.
    """
    if n <= 0:
        return 0.0, 0.0

    def chi2(p: float, df: int) -> float:
        if df == 0:
            return 0.0
        z = 1.959963984540054 * (1 if p > 0.5 else -1)
        return df * (1 - 2 / (9 * df) + z * math.sqrt(2 / (9 * df))) ** 3

    lo = chi2(0.025, 2 * k) / 2 if k else 0.0
    hi = chi2(0.975, 2 * (k + 1)) / 2
    return lo / n, hi / n


def _consensus_records(path: Path, min_depth: int, max_records: int):
    """Yield (sequence, depth) from a consensus FASTQ, keeping only groups of at least min_depth.

    The depth is the `cD:i:` tag assemble writes, which is the true read count of the molecule --
    not the number of reads that entered the consensus, which is capped.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    kept = 0
    with opener(path, "rt") as fh:
        while kept < max_records:
            head = fh.readline()
            if not head:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            fh.readline()
            depth = 0
            for field in head.split():
                if field.startswith("cD:i:"):
                    depth = int(field[5:])
                    break
            if depth >= min_depth:
                kept += 1
                yield seq, depth


def estimate_pre_amp_error(
    consensus_fastq: str | Path,
    min_depth: int = 20,
    max_minor: float = 0.01,
    max_divergence: float = 0.05,
    min_molecules: int = 100,
    min_bases: int = 100_000,
    max_records: int = 200_000,
) -> dict:
    """Fit the pre-amplification floor from deep consensuses, the way X2 measured it.

    A consensus over many reads has suppressed sequencing error to nothing, so what it still gets
    wrong is what was already in the molecule before the first amplification cycle. Compare deep
    consensuses against the library's own modal sequence and the residual IS the floor.

    Never: this needs a clonal template, and it refuses rather than guessing when it does not have
    one. On a diverse library every position is polymorphic, "disagrees with the modal base" means
    "is a different molecule", and the number that comes out is the library's diversity rather
    than its chemistry. Three exclusions, each of which showed up in the data (docs/nulls.rst,
    docs/quality_floor.rst):

      * a position whose minor allele reaches `max_minor` of molecules is real variation
      * a consensus disagreeing at more than `max_divergence` of positions is a different
        template -- an off-target product, not an erroneous copy
      * a tie (`N`) is not a call and does not vote

    Never: it is NOT a least-squares fit of `p_floor + a/c`. That model is wrong for a majority
    vote -- the sequencing residual falls geometrically in c, not as 1/c -- and on simulated data
    with a known floor it returned a negative probability.
    """
    path = Path(consensus_fastq)
    out: dict = {
        "min_depth": min_depth, "max_minor": max_minor, "max_divergence": max_divergence,
        "molecules_deep": 0, "molecules_scored": 0, "length": 0, "positions_scored": 0,
        "positions_polymorphic": 0, "mismatches": 0, "bases": 0,
        "rate": 0.0, "ci_lo": 0.0, "ci_hi": 0.0, "floor": 0.0, "ok": False, "reason": "",
    }
    seqs = [s for s, _ in _consensus_records(path, min_depth, max_records)]
    out["molecules_deep"] = len(seqs)
    if len(seqs) < min_molecules:
        out["reason"] = (
            f"only {len(seqs):,} molecules carry {min_depth} reads or more, against the "
            f"{min_molecules:,} it takes to see a floor. The floor is where the consensus curve "
            f"flattens, and a shallow library never gets there"
        )
        return out

    # One length, because a position only means anything if it is the same base in every molecule.
    # The modal length is the amplicon; anything else is a different product or a --contig
    # fragment, and mixing them would align position 40 of one against position 40 of another.
    lengths = collections.Counter(len(s) for s in seqs)
    length, at_length = lengths.most_common(1)[0]
    seqs = [s for s in seqs if len(s) == length]
    out["length"] = length
    if at_length < min_molecules:
        out["reason"] = (
            f"the modal consensus length {length} holds only {at_length:,} of "
            f"{out['molecules_deep']:,} deep molecules, so there is no common coordinate to score"
        )
        return out

    modal, minor = [], []
    for j in range(length):
        counts = collections.Counter(s[j] for s in seqs if s[j] != "N")
        total = sum(counts.values())
        base, n = counts.most_common(1)[0] if total else ("N", 0)
        modal.append(base)
        minor.append(0.0 if not total else 1.0 - n / total)
    scored = [j for j in range(length) if minor[j] < max_minor and modal[j] != "N"]
    out["positions_scored"] = len(scored)
    out["positions_polymorphic"] = length - len(scored)
    if len(scored) < length / 2:
        out["reason"] = (
            f"{out['positions_polymorphic']} of {length} positions carry a minor allele above "
            f"{max_minor:.1%}, so this library is not clonal enough to tell its own variation "
            f"from pre-amplification error. Name the chemistry instead"
        )
        return out

    mismatches = bases = 0
    for s in seqs:
        mm = sum(1 for j in scored if s[j] != modal[j] and s[j] != "N")
        called = sum(1 for j in scored if s[j] != "N")
        if called and mm / called > max_divergence:
            continue  # a different template, not an erroneous copy
        out["molecules_scored"] += 1
        mismatches += mm
        bases += called
    out["mismatches"], out["bases"] = mismatches, bases
    if out["molecules_scored"] < min_molecules or bases < min_bases:
        out["reason"] = (
            f"{out['molecules_scored']:,} molecules and {bases:,} scored bases survive the "
            f"exclusions, against {min_molecules:,} and {min_bases:,} needed. At a floor of 1e-4 "
            f"that is fewer than ten expected events, which is an interval and not a measurement"
        )
        return out

    out["rate"] = mismatches / bases
    out["ci_lo"], out["ci_hi"] = poisson_ci(mismatches, bases)
    # Zero observed errors is a bound, not a rate: take the upper end, because erring high costs a
    # few Phred points and erring low prints a confidence the data cannot support.
    floor = out["rate"] if mismatches else out["ci_hi"]
    if not 1e-7 <= floor <= 1e-2:
        out["reason"] = (
            f"the fitted floor {floor:.2e} is outside 1e-7..1e-2, which no chemistry is. Something "
            f"other than pre-amplification error is in this comparison"
        )
        return out
    out["floor"], out["ok"] = floor, True
    return out


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
    rt_floor: float | str = _core.RT_FLOOR,
    linkage_threshold: float = _core.LINKAGE_THRESHOLD,
    contig: bool = False,
    fast: bool = False,
    mate2: str | Path = "",
    merge_mates: bool = False,
    min_reads: int = 1,
    gzip_level: int = _core.GZIP_LEVEL,
    threads: int = 0,
    limit_reads: int = 0,
    limit_umis: int = 0,
) -> dict:
    """Assemble consensuses from a checkout-tagged FASTQ into `out_dir`.

    `reads` may also be the `.mig` buckets `checkout --mig` wrote -- one bucket file, a directory
    holding them, or a glob -- in which case the partition pass does not run at all: the reads are
    already range-partitioned on the key this stage groups by, which is the only thing that pass
    builds. One sample's buckets, never two: a UMI repeats across samples by design.

    `mate2` is the other mate of `reads` (checkout's `<sample>_R2.fq.gz`), matched by position;
    with `.mig` buckets the pair is already in the record and `merge_mates` alone turns it on.
    Either way mate 2 is reverse-complemented and PLACED against mate 1, so a pair whose mates
    overlap becomes one consensus spanning the insert and a pair whose mates do not stays two.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    estimate = None
    if isinstance(rt_floor, str):
        # `auto` costs a second assembly, deliberately. The consensus SEQUENCES do not depend on
        # the floor -- only the emitted quality does -- so the probe run is an ordinary assembly
        # whose output is read back and thrown away, and the real run then emits with the fitted
        # cap. The alternative was rewriting the qualities of a gzipped FASTQ in place, which
        # trades one clean pass for a new way to produce a broken file.
        probe = out / ".pre_amp_probe"
        probe_summary = run(
            reads, probe, sample_id=sample_id, rt_floor=_core.RT_FLOOR,
            linkage_threshold=linkage_threshold, contig=contig, fast=fast, mate2=mate2,
            merge_mates=merge_mates, min_reads=min_reads, gzip_level=gzip_level, threads=threads,
            limit_reads=limit_reads, limit_umis=limit_umis,
        )
        suffix = ".consensus.fq.gz" if gzip_level > 0 else ".consensus.fq"
        estimate = estimate_pre_amp_error(probe / (probe_summary["sample_id"] + suffix))
        shutil.rmtree(probe, ignore_errors=True)
        # Never falls back silently: the class is named in the report and written to the TSV next
        # to the numbers that refused it.
        rt_floor = estimate["floor"] if estimate["ok"] else RT_FLOORS["rt"]
        estimate["fallback"] = None if estimate["ok"] else "rt"
    mig_inputs = _mig_buckets(reads)
    summary = _core.assemble(
        "" if mig_inputs else str(reads), str(out), sample_id, rt_floor, linkage_threshold,
        contig, fast, min_reads, gzip_level, 0, threads, limit_reads, limit_umis, mig_inputs,
        str(mate2), merge_mates or bool(mate2),
    )
    summary["input"] = str(reads)
    summary["mig_input"] = mig_inputs
    summary["rt_floor"] = rt_floor
    summary["linkage_threshold"] = linkage_threshold
    if estimate is not None:
        summary["pre_amp_estimate"] = estimate
        # Never: a model-derived number is written next to what checks it. Every input the fit
        # rested on is in this row, so a floor that looks wrong can be argued with.
        with open(out / "assemble.pre_amp_error.tsv", "w") as fh:
            fh.write(
                "sample_id\tused\trate\tci_lo\tci_hi\tmismatches\tbases\tmolecules_scored\t"
                "molecules_deep\tmin_depth\tlength\tpositions_scored\tpositions_polymorphic\t"
                "fallback\treason\n"
            )
            e = estimate
            fh.write(
                f"{summary['sample_id']}\t{rt_floor:.3e}\t{e['rate']:.3e}\t{e['ci_lo']:.3e}\t"
                f"{e['ci_hi']:.3e}\t{e['mismatches']}\t{e['bases']}\t{e['molecules_scored']}\t"
                f"{e['molecules_deep']}\t{e['min_depth']}\t{e['length']}\t{e['positions_scored']}\t"
                f"{e['positions_polymorphic']}\t{e['fallback'] or ''}\t{e['reason']}\n"
            )

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
    if s.get("pre_amp_estimate"):
        e = s["pre_amp_estimate"]
        if e["ok"]:
            lines += [
                f"pre-amp floor fitted  {e['rate']:.2e} [{e['ci_lo']:.2e}, {e['ci_hi']:.2e}] from "
                f"{e['mismatches']:,} mismatches over {e['bases']:,} bases of "
                f"{e['molecules_scored']:,} molecules at {e['min_depth']}+ reads, "
                f"{e['positions_polymorphic']} of {e['length']} positions excluded as real "
                f"variation",
            ]
        else:
            lines += [
                f"pre-amp floor NOT fitted, using the '{e['fallback']}' class "
                f"({RT_FLOORS[e['fallback']]:.0e}): {e['reason']}",
            ]
    if s["fast_mode"]:
        lines += [
            f"modal sequence carried {s['mean_support']:.1%} of the reads in the groups that "
            f"emitted one -- counting mode, so no per-base correction and no sub-clustering",
        ]
    if s.get("merge_mates"):
        merged = s["groups"] - s["groups_fragmented"]
        lines += [
            "",
            f"mates merged  {merged:,} of {s['groups']:,} molecules came back as ONE contig "
            f"({_pct(merged, max(s['groups'], 1))}); the rest are pairs whose mates do not "
            f"overlap and are emitted as separate contigs rather than bridged",
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
