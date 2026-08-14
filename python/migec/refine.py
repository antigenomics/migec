"""The refine stage: correct barcode errors, then hand the reads on.

This is the stage that decides how many molecules there were. All of the work is in C++
(`_core.refine`); this module calls it once and turns the summary into tables and a report.
"""

from __future__ import annotations

import json
from pathlib import Path

from migec import _core
from migec.buckets import mig_buckets
from migec.checkout import _bytes, _dur, _pct


def run(
    reads: str | Path,
    out_dir: str | Path,
    sample_id: str = "",
    use_quality: bool = True,
    use_payload: bool = True,
    payload_width: int = 32,
    min_posterior: float = 0.95,
    expect_cells: int = 3000,
    cell_whitelist: str | Path = "",
    min_whitelist_posterior: float = 0.975,
    target_fdr: float = 0.05,
    gzip_level: int = _core.GZIP_LEVEL,
    threads: int = 0,
    limit_reads: int = 0,
    limit_umis: int = 0,
) -> dict:
    """Correct the barcodes in a checkout-tagged FASTQ, writing into `out_dir`.

    `reads` may also be the `.mig` buckets `checkout --mig` wrote -- one bucket file, or a
    directory holding one sample's. refine then reads those and writes `.mig` buckets back,
    re-partitioned on the CORRECTED barcode, which `assemble` takes directly. The audit trail is
    `<sample>.barcodes.tsv`: a `.mig` record has no room for the pre-correction barcode the way a
    FASTQ comment has `OX:Z:`, and the table carries it one row per barcode instead of once per
    read.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    mig_inputs = mig_buckets(reads)
    summary = _core.refine(
        "" if mig_inputs else str(reads), str(out), sample_id, use_quality, use_payload,
        payload_width, min_posterior, expect_cells, str(cell_whitelist), min_whitelist_posterior,
        target_fdr, gzip_level, threads, limit_reads, limit_umis, mig_inputs,
    )
    summary["input"] = str(reads)
    summary["mig_input"] = mig_inputs
    summary["min_posterior"] = min_posterior
    summary["target_fdr"] = target_fdr

    with open(out / "refine.coverage.tsv", "w") as fh:
        fh.write("sample_id\tmin_reads\tmax_reads\tmolecules\n")
        for b in summary["coverage"]:
            fh.write(
                f"{summary['sample_id']}\t{b['min_reads']}\t{b['max_reads']}\t{b['molecules']}\n"
            )
    (out / "refine.json").write_text(json.dumps(summary, indent=2, default=str))
    return summary


def format_report(summary: dict) -> str:
    s = summary
    lines = [
        f"reads       {s['reads']:,}",
        f"barcodes    {s['barcodes']:,} distinct",
        f"  merged    {s['merged']:,} ({_pct(s['merged'], max(s['barcodes'], 1))}) into a parent, "
        f"{s['merged_reads']:,} reads moved",
    ]
    if s["merged_by_payload"]:
        lines.append(
            f"    of which {s['merged_by_payload']:,} the count ratio alone would have refused "
            f"-- the reads agreed on the molecule"
        )
    lines += [
        "",
        f"molecules   {s['molecules']:,} after correction"
        + (
            f", {s['molecules_corrected']:,.0f} estimated after collisions"
            if not s["saturated"]
            else " (collision correction declined: the barcode space is saturated)"
        ),
        "",
    ]
    w = s["whitelist"]
    if w["barcodes"]:
        lines += [
            f"whitelist   {w['exact']:,} of {w['barcodes']:,} cell barcodes were on the list",
            f"            {w['corrected']:,} snapped to it ({w['reads_corrected']:,} reads), "
            f"{w['off_list']:,} left off it",
            f"            {w['far']:,} are further than one substitution from every entry -- "
            f"which is what measures the off-list prior, {w['background_prior']:.2e} per barcode",
            "",
        ]
    if s["cell_length"]:
        lines += [
            f"cells       {s['cells_called']:,} called of {s['cells_observed']:,} barcodes seen, "
            f"at >= {s['cell_threshold']:,} molecules (OrdMag)",
            f"            {s['molecules_in_called']:,} molecules in called cells "
            f"({_pct(s['molecules_in_called'], max(s['molecules'], 1))} of all)",
            f"            the curve breaks at rank {s['knee_rank']:,} "
            f"({s['knee_molecules']:,} molecules) -- the knee, for comparison",
            "",
        ]
    lines += [
        f"barcode error   {s['estimated_error']:.2e} per base, estimated from the "
        f"distance-1 excess",
        # The same rate off the children instead of off the neighbourhood, because the two fail in
        # different directions and quoting one alone hides which. Read at depth, since a child
        # whose parent was never sequenced cannot be counted and that is most of them at 1 read.
        f"                {s['error_at_depth']:.2e} from the reads in the children of molecules "
        f"seen >= {s['error_depth']} times"
        + (f" -- Q{s['error_phred']:.0f}" if s["error_phred"] else ""),
        f"                {s['error_from_children']:.2e} over all depths, which is a LOWER bound: "
        f"a child with no sequenced parent cannot be found",
        f"residual        {s['suspected_residual']:,} molecules still look like children of a "
        f"neighbour -- by count, or by their reads agreeing on the molecule",
        f"                {s['residual_fdr_at_one']:.2%} of 1-read molecules; "
        + (
            f"at >= {s['mig_size_threshold']:,} reads the rate is within the "
            f"{s['target_fdr']:.0%} target"
            if s["mig_size_threshold"]
            else f"no MIG size reaches the {s['target_fdr']:.0%} target"
        )
        + ". REPORTED, not applied -- a molecule seen three times is still a molecule",
        f"clonality       {s['payload_clonality']:.4f} of random barcode pairs carry the same "
        f"payload anyway",
        f"                -- payload agreement is worth about "
        f"{_worth(s['payload_clonality'])} here",
        "",
        f"{_dur(s['wall_seconds'])} = {_dur(s['table_seconds'])} table + "
        f"{_dur(s['correct_seconds'])} correction + {_dur(s['rewrite_seconds'])} rewrite, "
        f"on {s['threads']} threads",
        f"peak RSS {_bytes(s['peak_rss_bytes'])} of which the barcode table "
        f"{_bytes(s['table_bytes'])}",
        "",
        f"{'MIG size':>12}{'molecules':>12}{'share':>9}",
    ]
    total = sum(b["molecules"] for b in s["coverage"]) or 1
    for b in s["coverage"]:
        if not b["molecules"]:
            continue
        label = (
            f"{b['min_reads']}"
            if b["min_reads"] == b["max_reads"]
            else f"{b['min_reads']}-{b['max_reads']}"
        )
        lines.append(f"{label:>12}{b['molecules']:>12,}{100 * b['molecules'] / total:>8.1f}%")

    warnings = []
    if s.get("limited"):
        warnings.append(
            "the intake was limited (--limit-read / --limit-umi), so every number here describes "
            "the FIRST reads of the file rather than the library. Nothing measured under a limit "
            "-- error rate, occupancy, molecule count -- transfers to the whole run"
        )
    # OrdMag is a rule, not a measurement, and the knee is what the data says on its own. When
    # they disagree by more than a factor of three the rule is being applied to a curve it does
    # not describe -- an over-loaded run, or ambient RNA, or simply the wrong --expect-cells.
    if s["cell_length"] and s["cells_called"] and s["knee_rank"]:
        ratio = max(s["cells_called"], s["knee_rank"]) / max(
            min(s["cells_called"], s["knee_rank"]), 1
        )
        if ratio > 3:
            warnings.append(
                f"OrdMag calls {s['cells_called']:,} cells but the curve breaks at rank "
                f"{s['knee_rank']:,}, a factor of {ratio:.1f}. One of them is describing a "
                f"different library -- check --expect-cells and the rank plot before trusting "
                f"either"
            )
    singletons = next((b["molecules"] for b in s["coverage"] if b["min_reads"] == 1), 0)
    # The count ratio is the evidence a deep library has and a shallow one does not. Below ~3
    # reads/UMI most of what correction can do is done by the quality and the payload, and most
    # barcode errors are not fixable at all because the parent was never sequenced.
    if s["reads"] and s["barcodes"] and s["reads"] / s["barcodes"] < 3.0:
        warnings.append(
            f"{s['reads'] / s['barcodes']:.2f} reads per barcode. The count ratio carries almost "
            f"no evidence here, and a barcode error whose parent was never sequenced cannot be "
            f"corrected at all -- at ~1 read per barcode that is ~80% of them. What is corrected "
            f"is corrected conservatively; the rest inflates the molecule count rather than "
            f"corrupting a sequence"
        )
    if s["saturated"]:
        warnings.append(
            "the barcode space is saturated, so both the collision correction and the barcode "
            "error estimate are declined rather than guessed. A longer UMI is the fix"
        )
    if singletons > 0.5 * total:
        warnings.append(
            f"{_pct(singletons, total)} of molecules are still singletons after correction"
        )
    if warnings:
        lines.append("")
        lines += [f"warning: {w}" for w in warnings]
    return "\n".join(lines)


def _worth(clonality: float) -> str:
    """How much a matching payload is worth, in words rather than in a log."""
    if clonality <= 0.0:
        return "nothing -- payload evidence is off"
    if clonality > 0.5:
        return "nothing: this library is clonal, so two unrelated barcodes match anyway"
    return f"{1 / clonality:.0f}x odds towards the same molecule"
