"""The subsample stage: a smaller library that is still a library.

Sampling READS is trivial and every toolkit does it. Sampling MOLECULES is not, because a
molecule's reads are scattered through the file and nothing in a raw FASTQ says which they are --
which is why this exists. It is what equalises molecule counts across a cohort sequenced to
different depths, what a saturation curve is computed over, and what makes a test fixture behave
like the run it came from.

Never: Never a fraction of the reads. At four reads per molecule, ten thousand random reads give ten
thousand molecules seen once each -- the MIG size distribution is gone and every consensus is a
single read, so the smaller library is not a smaller version of the same library.
"""

from __future__ import annotations

from pathlib import Path

from migec import _core, bam
from migec.checkout import _bytes, _dur, _pct


def run(
    reads: str | Path,
    output: str | Path,
    keep_percent: float = 1.0,
    by_cell: bool = True,
    gzip_level: int = _core.GZIP_LEVEL,
) -> dict:
    """Keep all the reads of `keep_percent` of the barcodes.

    `reads` may be a BAM, SAM or CRAM carrying `RX` (see `migec.bam`); a paired one is sampled on
    mate 1.
    """
    per_10k = round(keep_percent * 100)
    if not 1 <= per_10k <= 10000:
        raise ValueError(
            f"--keep {keep_percent} is {per_10k} ten-thousandths; it must be in 0.01..100"
        )
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    if bam.is_alignment(reads):
        with bam.as_fastq(reads, Path(output).parent) as (mate1, mate2):
            summary = run(mate1, output, keep_percent, by_cell, gzip_level)
        summary["input"] = str(reads) + ("#R1" if mate2 is not None else "")
        return summary
    summary = _core.subsample(str(reads), str(output), per_10k, by_cell, gzip_level)
    summary["input"] = str(reads)
    summary["output"] = str(output)
    summary["keep_percent"] = keep_percent
    return summary


def format_report(summary: dict) -> str:
    s = summary
    barcodes = max(s["barcodes"], 1)
    lines = [
        f"read  {s['reads']:,}",
        f"kept  {s['reads_kept']:,} reads ({_pct(s['reads_kept'], max(s['reads'], 1))}) "
        f"in {s['barcodes']:,} barcodes",
        f"      {s['reads_kept'] / barcodes:.2f} reads per barcode on average, "
        f"median {s['reads_per_barcode_median']:,}, deepest {s['reads_per_barcode_max']:,} -- "
        f"the same distribution as the input, which is the point",
    ]
    if s["examples"]:
        shown = ", ".join(f"{bc} x{depth}" for bc, depth in s["examples"])
        lines.append(f"      {shown}")
        lines.append(
            "      (five kept barcodes in key order, not the first five seen: first-seen "
            "order is a sample of the deep MIGs and of nothing else)"
        )
    lines.append(
        f"{_dur(s['wall_seconds'])}, peak RSS {_bytes(s.get('peak_rss_bytes', 0))}"
    )
    if s["reads_without_umi"]:
        lines.append(
            f"warning: {s['reads_without_umi']:,} reads carried no RX tag and were dropped"
        )
    return "\n".join(lines)
