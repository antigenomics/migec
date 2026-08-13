"""The suggest stage: read the barcode layout off the reads.

A UMI cycle is one the synthesiser mixed: all four bases near 1/4, ~2 bits. A constant cycle is one
base at ~100%, ~0 bits. Everything else is payload. Segmenting the per-cycle composition on that
gives a pattern that can be pasted into a barcode table, which beats trusting a protocol
description written for the bench rather than for the file.
"""

from __future__ import annotations

import json
from pathlib import Path

from migec import _core


def run(
    reads: str | Path,
    out_dir: str | Path | None = None,
    cycles: int = 60,
    max_reads: int = 200_000,
    umi_deviation: float = 0.18,
) -> dict:
    """Profile `reads` and suggest a barcode pattern. Writes TSVs if `out_dir` is given."""
    summary = _core.suggest(str(reads), cycles, max_reads, umi_deviation)
    summary["input"] = str(reads)
    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        with open(out / "suggest.cycles.tsv", "w") as fh:
            fh.write(
                "cycle\tA\tC\tG\tT\tentropy_bits\tcollision\tconsensus\tconsensus_fraction\t"
                "deviation_from_uniform\tmean_phred\n"
            )
            for c in summary["cycles"]:
                fh.write(
                    f"{c['cycle']}\t{c['A']:.6f}\t{c['C']:.6f}\t{c['G']:.6f}\t{c['T']:.6f}\t"
                    f"{c['entropy']:.6f}\t{c['collision']:.6f}\t{c['consensus']}\t"
                    f"{c['consensus_fraction']:.6f}\t{c['deviation']:.6f}\t{c['mean_phred']:.2f}\n"
                )
        with open(out / "suggest.segments.tsv", "w") as fh:
            fh.write("kind\tbegin\tend\tlength\tconsensus\tmean_deviation\n")
            for s in summary["segments"]:
                fh.write(
                    f"{s['kind']}\t{s['begin']}\t{s['end']}\t{s['length']}\t{s['consensus']}\t"
                    f"{s['mean_deviation']:.6f}\n"
                )
        (out / "suggest.json").write_text(json.dumps(summary, indent=2, default=str))
    return summary


def format_report(summary: dict) -> str:
    lines = [
        f"profiled {summary['reads']:,} reads, {summary['read_length']} nt",
        "",
        f"{'cycle':>6}{'A':>7}{'C':>7}{'G':>7}{'T':>7}{'1/4 dev':>9}{'Q':>6}  layout",
    ]
    kind_of = {}
    for s in summary["segments"]:
        for j in range(s["begin"], s["end"]):
            kind_of[j] = s["kind"]
    glyph = {"umi": "N  UMI", "constant": "|  constant", "variable": ".  variable"}
    shown = 0
    for c in summary["cycles"]:
        k = kind_of.get(c["cycle"], "variable")
        # Constant and payload runs are long and uninformative past the first few cycles.
        first_of_run = kind_of.get(c["cycle"] - 1) != k
        if k == "umi" or first_of_run or shown < 2:
            lines.append(
                f"{c['cycle']:>6}{c['A']:>7.3f}{c['C']:>7.3f}{c['G']:>7.3f}{c['T']:>7.3f}"
                f"{c['deviation']:>9.3f}{c['mean_phred']:>6.0f}  {glyph[k]}"
            )
            shown = 0 if first_of_run else shown + 1
        elif shown == 2:
            lines.append(f"{'...':>6}")
            shown += 1

    lines += ["", "segments:"]
    for s in summary["segments"]:
        extra = f"  {s['consensus']}" if s["kind"] == "constant" else ""
        lines.append(
            f"  {s['begin']:>3}-{s['end'] - 1:<3} {s['kind']:<9} {s['length']:>3} nt"
            f"  (mean 1/4 deviation {s['mean_deviation']:.3f}){extra}"
        )

    lines += [
        "",
        f"UMI      {summary['umi_length']} nt   nominal space 4^{summary['umi_length']} = "
        f"{4 ** summary['umi_length']:,}",
        f"anchor   {summary['anchor_length']} nt of constant sequence to score against",
        "",
        "pattern  " + (summary["pattern"] or "(none)"),
    ]
    if summary["umi_length"]:
        # A real tab, because that is what a barcode table needs and this line is meant to be
        # copied out of the terminal into one.
        lines.append(f"\npaste into a barcode table as:\n  S1\t{summary['pattern']}")
    if summary["note"]:
        lines.append(f"\nwarning: {summary['note']}")
    return "\n".join(lines)
