"""The suggest stage: read the barcode layout off the reads.

A UMI cycle is one the synthesiser mixed: all four bases near 1/4, ~2 bits. A constant cycle is one
base at ~100%, ~0 bits. Everything else is payload. Segmenting the per-cycle composition on that
gives a pattern that can be pasted into a barcode table, which beats trusting a protocol
description written for the bench rather than for the file.
"""

from __future__ import annotations

import json
from pathlib import Path

from migec import _core, bam


def run(
    reads: str | Path,
    out_dir: str | Path | None = None,
    cycles: int = 60,
    max_reads: int = 200_000,
    umi_deviation: float = 0.18,
) -> dict:
    """Profile `reads` and suggest a barcode pattern. Writes TSVs if `out_dir` is given.

    `reads` may be a BAM, SAM or CRAM. No `RX` is required -- the whole point of this stage is that
    nothing is known about the layout yet -- and a paired file is profiled on mate 1, the same
    convention as a paired FASTQ, which this stage also takes one mate at a time.
    """
    if bam.is_alignment(reads):
        with bam.as_fastq(reads, out_dir or Path.cwd(), need_rx=False) as (mate1, mate2):
            summary = run(mate1, out_dir, cycles, max_reads, umi_deviation)
        summary["input"] = str(reads) + ("#R1" if mate2 is not None else "")
        if out_dir is not None:
            (Path(out_dir) / "suggest.json").write_text(json.dumps(summary, indent=2, default=str))
        return summary
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
        with open(out / "suggest.kmers.tsv", "w") as fh:
            fh.write("kmer\tcount\texpected\tratio\tmean_position\tread_fraction\n")
            for k in summary["kmers"]:
                fh.write(
                    f"{k['kmer']}\t{k['count']}\t{k['expected']:.1f}\t{k['ratio']:.3f}\t"
                    f"{k['mean_position']:.1f}\t{k['read_fraction']:.6f}\n"
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
    # Overrepresented k-mers: what synthetic sequence is still in these reads. On raw input that
    # is the primer the pattern is about to be built from; run on trimmed or consensus output it
    # is whatever the trim failed to remove, which is the only way to find out.
    strong = [k for k in summary["kmers"] if k["ratio"] >= 5.0 and k["read_fraction"] >= 0.01]
    if strong:
        lines += [
            "",
            f"{'kmer':<10}{'count':>10}{'obs/exp':>10}{'reads':>9}{'mean pos':>10}",
        ]
        for k in strong[:10]:
            lines.append(
                f"{k['kmer']:<10}{k['count']:>10,}{k['ratio']:>10.1f}"
                f"{k['read_fraction']:>8.1%}{k['mean_position']:>10.1f}"
            )
        joined = _stitch(strong)
        if joined:
            lines.append(f"\noverlapping into: {joined}")

    if summary["note"]:
        lines.append(f"\nwarning: {summary['note']}")
    return "\n".join(lines)


def _stitch(hits: list[dict], min_count_ratio: float = 0.5) -> str:
    """Greedily join overlapping k-mers back into the sequence they came from, both ways.

    Eight bases name a primer but do not print one. Real synthetic sequence shows up as a run of
    k-mers each shifted one base from the last, so following that chain out of the strongest hit
    recovers the adapter itself, which is what the reader needs in order to act on it.

    The chain stops where the count drops: a k-mer straddling the end of the primer is carried by
    whatever follows it, so it is seen a quarter as often per base of overhang. Without that test
    the walk keeps going one base at a time into the payload and prints an adapter that is partly
    invented.
    """
    counts = {k["kmer"]: k["count"] for k in hits}
    seed = hits[0]["kmer"]
    floor = min_count_ratio * counts[seed]
    k = len(seed)

    def extend(seq: str, table: dict[str, int]) -> str:
        seen = {seq[-k:]}
        for _ in range(120):  # bounded: a repeat would otherwise walk forever
            nxt = max(
                (m for m in table
                 if m.startswith(seq[-(k - 1):]) and table[m] >= floor and m not in seen),
                key=table.get,
                default=None,
            )
            if nxt is None:
                break
            seen.add(nxt)
            seq += nxt[-1]
        return seq

    out = extend(seed, counts)
    reversed_counts = {m[::-1]: c for m, c in counts.items()}
    out = extend(out[::-1], reversed_counts)[::-1]
    return out if len(out) > k else ""
