"""The assemble stage: reads sharing a UMI are reads of one molecule, so collapse them.

All of the work is in C++ (`_core.assemble`); this module calls it once and turns the returned
summary into tables and a report. The per-molecule table is written by the C++ side because it is
one line per molecule and there are hundreds of millions of them.
"""

from __future__ import annotations

import json
from pathlib import Path

from migec import _core
from migec.checkout import _bytes, _dur, _pct


def run(
    reads: str | Path,
    out_dir: str | Path,
    sample_id: str = "",
    rt_floor: float = 1e-4,
    linkage_threshold: float = 9.61,
    contig: bool = False,
    min_reads: int = 1,
    gzip_level: int = 6,
) -> dict:
    """Assemble consensuses from a checkout-tagged FASTQ into `out_dir`."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary = _core.assemble(
        str(reads), str(out), sample_id, rt_floor, linkage_threshold, contig, min_reads,
        gzip_level
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
        f"{s['buckets']} bucket(s) + {_dur(s['wall_seconds'] - s['partition_seconds'])} consensus",
        f"peak RSS {_bytes(s['peak_rss_bytes'])} -- one bucket is resident at a time",
        "",
        f"mean emitted quality  Q{s['mean_quality']:.1f}  "
        f"(capped at Q{s['quality_cap']:.0f} by the RT floor of {s['rt_floor']:.1e})",
        f"mean consensus error  {s['mean_consensus_error']:.2e} before the floor is added",
        "",
        f"{'MIG size':>12}{'groups':>12}{'share':>9}",
    ]
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
            f"molecule. X3 measured a 1.26% rate on a real library at this threshold "
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
