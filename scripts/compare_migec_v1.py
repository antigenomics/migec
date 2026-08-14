#!/usr/bin/env python3
# 2026-08-14
#
# migec 2 against MIGEC 1.2.9, the Groovy implementation this repo replaced.
#
# The rewrite's whole claim is that the algorithms are the same specification and the
# implementation is not, so the comparison that matters is: does it find the same molecules, and
# how much faster. Both are measured here on one simulated library with a known truth --
# `n_molecules` templates, so a molecule count that is right is a count near it and a consensus
# that is right is one that matches `truth_consensus.fa`.
#
# Note: the two pipelines are not the same shape and the mapping is stated rather than assumed.
# v1 is Checkout -> Assemble, with barcode-error correction folded into Assemble as
# `--filter-collisions` (a count-ratio rule). v2 is checkout -> refine -> assemble, with correction
# its own stage and its own evidence. `--filter-collisions` is therefore ON for v1, or v1 is being
# asked to do less work and the wall-clock ratio flatters us.
#
# Usage:
#     python scripts/compare_migec_v1.py --out /tmp/v1 --jar migec-1.2.9.jar --molecules 20000
#
# Get the jar: gh release download 1.2.9 --repo antigenomics/migec -p 'migec-1.2.9.zip'

from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def _run(argv: list[str]) -> tuple[float, int]:
    """Wall clock and peak RSS of one child, from wait4's own rusage."""
    t0 = time.perf_counter()
    p = subprocess.Popen(argv, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    _, status, ru = os.wait4(p.pid, 0)
    if status:
        raise SystemExit(f"failed ({status}): {' '.join(argv)}")
    scale = 1 if sys.platform == "darwin" else 1024
    return time.perf_counter() - t0, ru.ru_maxrss * scale


def read_consensus(path: Path) -> list[str]:
    """Consensus sequences from a FASTQ, uppercased. Names differ between the two tools."""
    out = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                out.append(line.strip().upper())
    return out


def read_truth_consensus(path: Path) -> set[str]:
    seqs, cur = set(), []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.add("".join(cur).upper())
                cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.add("".join(cur).upper())
    return seqs


def score(consensus: list[str], truth: set[str]) -> dict:
    """How many consensuses are exactly a template, and how many templates were recovered.

    Exact match, not an alignment: a consensus that is one base off is a consensus that would put
    a false variant into a caller, which is the whole thing the stage exists to prevent.
    """
    exact = sum(1 for s in consensus if s in truth)
    return {
        "consensuses": len(consensus),
        "exact": exact,
        "precision": exact / len(consensus) if consensus else 0.0,
        "recall": len(set(consensus) & truth) / len(truth) if truth else 0.0,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--jar", type=Path, required=True, help="migec-1.2.9.jar")
    ap.add_argument("--molecules", type=int, default=20_000)
    ap.add_argument("--clones", type=int, default=200)
    ap.add_argument("--coverage", type=float, default=8.0)
    ap.add_argument("--umi-len", type=int, default=12)
    ap.add_argument("--umi-error", type=float, default=2e-3)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--java", default="java")
    ap.add_argument("--min-count", type=int, default=1,
                    help="minimal reads per MIG, applied to BOTH. Never: v1 defaults to 5 and "
                         "migec 2 to 1, so leaving each at its own default compares defaults "
                         "rather than algorithms -- v1's output file is named .t5. for exactly "
                         "this reason")
    ap.add_argument("--tsv", type=Path)
    args = ap.parse_args(argv)

    from tests.synthetic._sim import SimConfig, simulate

    args.out.mkdir(parents=True, exist_ok=True)
    cfg = SimConfig(n_molecules=args.molecules, n_clones=args.clones, coverage=args.coverage,
                    coverage_cv=0.4, umi_len=args.umi_len, umi_error=args.umi_error,
                    adapter=ADAPTER, seed=args.seed)
    sim = simulate(cfg, args.out / "sim")
    truth = read_truth_consensus(Path(sim["truth_consensus"]))
    print(f"# {sim['n_reads']:,} reads, {sim['n_molecules']:,} molecules, "
          f"{len(truth):,} distinct templates", file=sys.stderr)

    sheet = args.out / "bc.txt"
    # v1's barcode file is `sample<TAB>master<TAB>slave`; a dot is "no slave".
    sheet.write_text(f"S1\t{sim['pattern']}\t.\n")

    rows = []

    # ---------------------------------------------------------------- MIGEC 1.2.9
    v1 = args.out / "v1"
    (v1 / "co").mkdir(parents=True, exist_ok=True)
    (v1 / "asm").mkdir(parents=True, exist_ok=True)
    t, rss = _run([args.java, "-jar", str(args.jar), "Checkout", "-cute", str(sheet),
                   sim["reads"], ".", str(v1 / "co")])
    # Note: v1 names a single-end checkout `S1_R0.fastq.gz` and a paired one `S1_R1/R2` -- the
    # mate number, not a file index -- so the file is found rather than named.
    co_out = next(iter(sorted((v1 / "co").glob("S1_R*.fastq.gz"))), None)
    if co_out is None:
        raise SystemExit(f"MIGEC v1 Checkout assigned no reads into {v1 / 'co'}")
    t2, rss2 = _run([args.java, "-jar", str(args.jar), "Assemble", "-c", "--filter-collisions",
                     "-m", str(args.min_count), str(co_out), ".", str(v1 / "asm")])
    v1_out = next(iter(sorted((v1 / "asm").glob("*.fastq.gz"))), None)
    if v1_out is None:
        raise SystemExit(f"MIGEC v1 wrote no consensus into {v1 / 'asm'}")
    rows.append({"tool": "migec-1.2.9", "seconds": t + t2, "peak_rss_bytes": max(rss, rss2),
                 **score(read_consensus(v1_out), truth)})

    # ---------------------------------------------------------------- migec 2
    from migec.assemble import run as assemble_run
    from migec.checkout import run as checkout_run
    from migec.refine import run as refine_run

    v2 = args.out / "v2"
    t0 = time.perf_counter()
    checkout_run(sim["reads"], sheet, v2 / "co", threads=args.threads)
    refine_run(v2 / "co" / "S1.fq.gz", v2 / "ref", threads=args.threads)
    st = assemble_run(v2 / "ref" / "S1.fq.gz", v2 / "asm", threads=args.threads,
                      min_reads=args.min_count)
    seconds = time.perf_counter() - t0
    v2_out = next(iter(sorted((v2 / "asm").glob("*.fq.gz"))), None)
    if v2_out is None:
        raise SystemExit(f"migec 2 wrote no consensus into {v2 / 'asm'}")
    rows.append({"tool": "migec-2", "seconds": seconds,
                 "peak_rss_bytes": st.get("peak_rss_bytes", 0),
                 **score(read_consensus(v2_out), truth)})

    cols = ["tool", "consensuses", "exact", "precision", "recall", "seconds", "peak_rss_bytes"]
    lines = ["\t".join(cols)]
    for r in rows:
        lines.append("\t".join(
            str(r[c]) if c == "tool" else
            f"{r[c]:.0f}" if c in ("consensuses", "exact", "peak_rss_bytes") else f"{r[c]:.4f}"
            for c in cols))
    speedup = rows[0]["seconds"] / rows[1]["seconds"] if rows[1]["seconds"] else 0.0
    out = "\n".join(lines)
    print(out)
    print(f"# migec 2 is {speedup:.1f}x MIGEC 1.2.9's wall clock", file=sys.stderr)
    if args.tsv:
        args.tsv.write_text(out + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
