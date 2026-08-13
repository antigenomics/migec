#!/usr/bin/env python
"""Thread scaling for `migec checkout`, measured rather than quoted.

2026-08-13. The README carried one set of numbers and ROADMAP.md another, which is what happens
when a table is edited by hand after a run. This writes the table, so the figure and the prose
come from the same measurement:

    python scripts/benchmark_threads.py --reads 2000000 -o assets/
    migec plot assets/            # draws benchmark_threads.svg from the TSV it just wrote

Both columns are reported because they scale with different things. `match_seconds` covers the
parallel driver -- read a chunk, match it, compress it, append it -- and threads. The end-to-end
clock also covers the per-sample statistics, which run once, on one thread, over every DISTINCT
barcode, and do not.
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
from pathlib import Path

from migec.checkout import run

BARCODES = """\
S1\taaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
S2\taaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN
S3\taaGCCcagtggtatcaacgcagagtNNNNtNNNNtNNNN
S4\taaGGTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
"""
TAGS = {"S1": "ACT", "S2": "AGA", "S3": "GCC", "S4": "GGT"}
ADAPTER = "CAGTGGTATCAACGCAGAGT"


def corpus(path: Path, n_reads: int, reads_per_umi: int, payload: int) -> Path:
    """The same shape tests/benchmark uses: four patterns, a 12 nt split UMI, shallow coverage."""
    reads = path / "reads.fq.gz"
    if reads.exists():
        return reads
    rng = random.Random(0)
    with gzip.open(reads, "wt", compresslevel=1) as fh:
        i = 0
        while i < n_reads:
            # On the MOLECULE, not on the read counter. `i % 4` with four reads per molecule is
            # always 0, so every read goes to S1 and three of the four patterns never match --
            # which still measures the matcher, and measures nothing about demultiplexing.
            sample = list(TAGS)[(i // reads_per_umi) % 4]
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            body = "".join(rng.choice("ACGT") for _ in range(payload))
            for _ in range(reads_per_umi):
                seq = ("AA" + TAGS[sample] + ADAPTER + umi[:4] + "T" + umi[4:8] + "T" + umi[8:]
                       + body)
                fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
                i += 1
    (path / "barcodes.txt").write_text(BARCODES)
    return reads


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reads", type=int, default=2_000_000)
    ap.add_argument("--reads-per-umi", type=int, default=4)
    ap.add_argument("--payload", type=int, default=90)
    ap.add_argument("--threads", default="1,2,4,8,16")
    ap.add_argument("-o", "--out", type=Path, default=Path("assets"))
    ap.add_argument("--work", type=Path, default=Path("scratch/benchmark"))
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    reads = corpus(args.work, args.reads, args.reads_per_umi, args.payload)

    rows = []
    for t in [int(x) for x in args.threads.split(",")]:
        s = run(reads, args.work / "barcodes.txt", args.work / f"out{t}", threads=t)
        match = s["match_seconds"] or s["wall_seconds"]
        rows.append({
            "threads": s["threads"],
            "reads_per_second": s["total"] / s["wall_seconds"],
            "matching_reads_per_second": s["total"] / match,
            "peak_rss_mb": s["peak_rss_bytes"] / 1e6,
            "wall_seconds": s["wall_seconds"],
            "match_seconds": match,
        })
        print(f"{t:>3} threads  {rows[-1]['reads_per_second']:>10,.0f} reads/s end to end  "
              f"{rows[-1]['matching_reads_per_second']:>10,.0f} matching  "
              f"{rows[-1]['peak_rss_mb']:>7.0f} MB")

    tsv = args.out / "benchmark_threads.tsv"
    with open(tsv, "w") as fh:
        fh.write("threads\treads_per_second\tmatching_reads_per_second\tpeak_rss_mb\t"
                 "wall_seconds\tmatch_seconds\n")
        for r in rows:
            fh.write(f"{r['threads']}\t{r['reads_per_second']:.0f}\t"
                     f"{r['matching_reads_per_second']:.0f}\t{r['peak_rss_mb']:.1f}\t"
                     f"{r['wall_seconds']:.3f}\t{r['match_seconds']:.3f}\n")
    (args.out / "benchmark_threads.json").write_text(json.dumps(
        {"reads": args.reads, "reads_per_umi": args.reads_per_umi, "payload_nt": args.payload,
         "rows": rows}, indent=2))
    print(f"\nwrote {tsv}")


if __name__ == "__main__":
    main()
