"""The checkout stage: demultiplex, extract UMIs, trim, and write the QC tables.

Everything per-read happens in C++ (`_core.run_checkout`); this module reads the sample sheet,
calls it once, and turns the returned summary into tables. Plots are generated from those tables,
never from inside the C++ -- so a figure can always be redrawn from a committed TSV.
"""

from __future__ import annotations

import json
from pathlib import Path

from migec import _core
from migec.sheet import SampleRow, read_barcodes


def run(
    reads: str | Path,
    barcodes: str | Path,
    out_dir: str | Path,
    reads2: str | Path | None = None,
    trim: str = "pattern",
    min_umi_quality: int = 0,
    write_unmatched: bool = False,
    threads: int = 0,
) -> dict:
    """Demultiplex `reads` (and `reads2`, if paired) using `barcodes`, writing into `out_dir`."""
    rows: list[SampleRow] = read_barcodes(barcodes)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    summary = _core.run_checkout(
        str(reads),
        "" if reads2 is None else str(reads2),
        [r.sample_id for r in rows],
        [r.pattern for r in rows],
        str(out) + "/",
        trim,
        min_umi_quality,
        write_unmatched,
        threads,
    )
    summary["input"] = str(reads)
    summary["input2"] = "" if reads2 is None else str(reads2)
    summary["barcodes"] = str(barcodes)
    summary["patterns"] = {r.sample_id: r.pattern for r in rows}

    _write_tables(out, summary)
    (out / "checkout.json").write_text(json.dumps(summary, indent=2, default=str))
    return summary


def _write_tables(out: Path, summary: dict) -> None:
    # One row per sample: what came out, and whether it is worth assembling.
    with open(out / "checkout.summary.tsv", "w") as fh:
        fh.write(
            "sample_id\treads\tumis\tmean_reads_per_umi\treads_in_migs_ge5\tover_sequenced\t"
            "umi_length\teffective_length\teffective_space\ttotal_entropy\ttotal_information\t"
            "umi_error_rate\tumis_merged\treads_merged\tmolecules_observed\tmolecules_corrected\t"
            "saturated\n"
        )
        for s in summary["samples"]:
            fh.write(
                f"{s['sample_id']}\t{s['reads']}\t{s['umis']}\t{s['mean_reads_per_umi']:.4f}\t"
                f"{s['reads_in_migs_ge5']:.6f}\t{int(s['over_sequenced'])}\t{s['umi_length']}\t"
                f"{s['effective_length']:.4f}\t{s['effective_space']:.1f}\t"
                f"{s['total_entropy']:.4f}\t{s['total_information']:.4f}\t"
                f"{s['umi_error_rate']:.3e}\t{s['umis_merged']}\t{s['reads_merged']}\t"
                f"{s['molecules_observed']}\t{s['molecules_corrected']:.1f}\t"
                f"{int(s['saturated'])}\n"
            )

    # The coverage histogram, in MIGEC's power-of-two bins so published figures are comparable.
    with open(out / "checkout.coverage.tsv", "w") as fh:
        fh.write("sample_id\tmig_size\treads\tunits\n")
        for s in summary["samples"]:
            for b, (reads, units) in enumerate(zip(s["hist_reads"], s["hist_units"])):
                if reads or units:
                    fh.write(f"{s['sample_id']}\t{2**b}\t{reads}\t{units}\n")

    # Per-position base usage and information content -- the numbers a sequence logo draws.
    with open(out / "checkout.umi_composition.tsv", "w") as fh:
        fh.write("sample_id\tposition\tA\tC\tG\tT\tentropy_bits\tinformation_bits\tcollision\n")
        for s in summary["samples"]:
            for p in s["composition"]:
                fh.write(
                    f"{s['sample_id']}\t{p['position']}\t{p['A']:.6f}\t{p['C']:.6f}\t"
                    f"{p['G']:.6f}\t{p['T']:.6f}\t{p['entropy']:.6f}\t{p['information']:.6f}\t"
                    f"{p['collision']:.6f}\n"
                )


def format_report(summary: dict) -> str:
    """A short human-readable report, printed at the end of a run."""
    c = summary
    lines = [
        f"reads       {c['total']:,}",
        f"  assigned  {c['assigned']:,} ({_pct(c['assigned'], c['total'])})",
        f"  unmatched {c['unmatched']:,} ({_pct(c['unmatched'], c['total'])})",
        f"  ambiguous {c['ambiguous']:,} ({_pct(c['ambiguous'], c['total'])})",
    ]
    if c["bad_umi"]:
        lines.append(f"  bad UMI   {c['bad_umi']:,} ({_pct(c['bad_umi'], c['total'])})")
    if c["short_payload"]:
        lines.append(f"  too short {c['short_payload']:,} ({_pct(c['short_payload'], c['total'])})")
    if c.get("normalised"):
        lines.append(
            f"  flipped   {c['normalised']:,} ({_pct(c['normalised'], c['total'])}) "
            f"-- tag found on the other mate/strand, reads normalised"
        )
    lines.append("")
    lines.append(
        f"{_dur(c.get('wall_seconds', 0.0))} "
        f"({c.get('reads_per_second', 0.0):,.0f} reads/s on {c.get('threads', 1)} threads), "
        f"peak RSS {_bytes(c.get('peak_rss_bytes', 0))} "
        f"of which UMI counters {_bytes(c.get('umi_memory_bytes', 0))}"
    )
    lines.append("")
    lines.append(
        f"{'sample':<12}{'reads':>12}{'UMIs':>12}{'reads/UMI':>11}{'UMI len':>9}{'eff len':>9}"
    )
    for s in summary["samples"]:
        lines.append(
            f"{s['sample_id']:<12}{s['reads']:>12,}{s['umis']:>12,}"
            f"{s['mean_reads_per_umi']:>11.2f}{s['umi_length']:>9}{s['effective_length']:>9.2f}"
        )

    warnings = []
    # The UMI counters are the one allocation that grows with the library rather than with the
    # chunk size, so they are the thing that decides whether a run fits. Range partitioning is
    # what fixes it; until that lands (M2) the honest thing is to say so before the OOM.
    if c.get("umi_memory_bytes", 0) > 1 << 30:
        warnings.append(
            f"UMI counters hold {_bytes(c['umi_memory_bytes'])}. This grows with the number of "
            f"distinct UMIs and is not yet partitioned across buckets, so a much larger input may "
            f"not fit in memory"
        )
    if c["total"] and c["assigned"] / c["total"] < 0.5:
        warnings.append(
            "less than half of reads matched a pattern -- run `migec suggest` to check where the "
            "barcode actually is"
        )
    for s in summary["samples"]:
        if not s["over_sequenced"]:
            warnings.append(
                f"{s['sample_id']}: {s['mean_reads_per_umi']:.1f} reads/UMI. Consensus assembly "
                f"needs over-sequencing; below ~5 most molecules are seen once"
            )
        if s["umi_length"] and s["effective_length"] < 0.8 * s["umi_length"]:
            warnings.append(
                f"{s['sample_id']}: UMI is {s['umi_length']} nt but only "
                f"{s['effective_length']:.1f} nt of it is usable -- the base composition is skewed, "
                f"so collisions are more frequent than the length suggests"
            )
        if s["saturated"]:
            warnings.append(
                f"{s['sample_id']}: observed UMIs are a large fraction of the usable space; "
                f"molecule counts are biased low and correction is deliberately conservative"
            )
    if warnings:
        lines.append("")
        for w in warnings:
            lines.append(f"warning: {w}")
    return "\n".join(lines)


def _pct(n: int, total: int) -> str:
    return f"{100.0 * n / total:.1f}%" if total else "n/a"


def _bytes(n: int) -> str:
    for unit in ("B", "kB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} GB"


def _dur(s: float) -> str:
    if s < 60:
        return f"{s:.1f} s"
    return f"{int(s // 60)}m{s % 60:04.1f}s"
