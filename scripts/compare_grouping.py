#!/usr/bin/env python3
# 2026-08-14
#
# UMI grouping accuracy, wall clock and peak memory against the map-first tools: UMI-tools
# (Smith et al., Genome Res 2017) and fgbio GroupReadsByUmi (Fulcrum Genomics).
#
# All three answer the same question -- which reads came from the same original molecule -- so the
# comparison is a clustering comparison against the simulator's known truth, scored with the
# adjusted Rand index by `compare_calib.adjusted_rand`. Calib has its own driver
# (`compare_calib.py`) because it takes raw FASTQ; these two do not.
#
# Note: the two of them group on (POSITION, UMI) and so cannot run at all without an alignment.
# That is the difference the numbers are about rather than an inconvenience: on a repertoire
# library everything maps to the same few references and the position carries almost no key bits,
# which is exactly the regime `clones.fa` simulates. `docs/downstream.rst` has the argument; this
# script has the measurement.
#
# Usage:
#     python scripts/compare_grouping.py --out /tmp/cmp --molecules 20000 --coverage 5
#     python scripts/compare_grouping.py --out /tmp/cmp --tools migec,umi_tools   # skip fgbio
#
# Needs on PATH: minimap2, samtools, and whichever of `umi_tools` / `fgbio` is asked for. fgbio
# needs a JDK 17 or newer; point JAVA_HOME at one if the default is older.

from __future__ import annotations

import argparse
import gzip
import os
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from compare_calib import adjusted_rand, read_migec, read_truth  # noqa: E402

ADAPTER = "CAGTGGTATCAACGCAGAGT"


class Cost:
    """Wall clock and peak RSS of a subprocess tree, from wait4's own rusage.

    Never: `getrusage(RUSAGE_CHILDREN)` is a high-water mark over EVERY child this process has
    reaped, so a difference between two calls is not one tool's memory -- run the second tool
    first and it reads zero. wait4 reports the child that just exited.
    """

    def __init__(self) -> None:
        self.seconds = 0.0
        self.rss = 0

    def run(self, argv: list[str], stdout=None, env: dict | None = None) -> None:
        t0 = time.perf_counter()
        p = subprocess.Popen(argv, stdout=stdout, stderr=subprocess.DEVNULL,
                             env={**os.environ, **(env or {})})
        _, status, ru = os.wait4(p.pid, 0)
        self.seconds += time.perf_counter() - t0
        # macOS reports bytes, Linux kilobytes.
        scale = 1 if sys.platform == "darwin" else 1024
        self.rss = max(self.rss, ru.ru_maxrss * scale)
        if status:
            raise SystemExit(f"{argv[0]} failed with status {status}: {' '.join(argv)}")


def simulate(out: Path, molecules: int, clones: int, coverage: float, umi_len: int,
             umi_error: float, seed: int) -> dict:
    from tests.synthetic._sim import SimConfig, simulate as sim

    cfg = SimConfig(n_molecules=molecules, n_clones=clones, coverage=coverage, coverage_cv=0.4,
                    umi_len=umi_len, umi_error=umi_error, adapter=ADAPTER, seed=seed)
    return sim(cfg, out / "sim")


def run_migec(sim: dict, out: Path, threads: int) -> tuple[dict[str, str], Cost]:
    """checkout then refine. The group is the CORRECTED barcode, which is what refine is for."""
    from migec.checkout import run as checkout_run
    from migec.refine import run as refine_run

    sheet = out / "bc.txt"
    sheet.write_text(f"S1\t{sim['pattern']}\n")
    cost = Cost()
    t0 = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    checkout_run(sim["reads"], sheet, out / "co", threads=threads)
    stats = refine_run(out / "co" / "S1.fq.gz", out / "ref", threads=threads)
    cost.seconds = time.perf_counter() - t0
    scale = 1 if sys.platform == "darwin" else 1024
    cost.rss = max(stats["peak_rss_bytes"], before * scale)
    return read_migec([out / "ref" / "S1.fq.gz"]), cost


def _strip_umi_into_name(src: str, dst: Path, umi_len: int) -> None:
    """UMI-tools' own convention: the barcode moves into the read name after an underscore.

    This is `umi_tools extract` done inline, because extract compiles a pattern and we already
    know the layout -- and because doing it here keeps its cost out of the timing.
    """
    with gzip.open(src, "rt") as fh, gzip.open(dst, "wt") as o:
        name = ""
        for i, line in enumerate(fh):
            if i % 4 == 0:
                name = line[1:].split()[0]
            elif i % 4 == 1:
                o.write(f"@{name}_{line[:umi_len]}\n{line[umi_len:]}")
            elif i % 4 == 2:
                o.write("+\n")
            else:
                o.write(line[umi_len:])


def _align(reads: Path, ref: str, bam: Path, cost: Cost, threads: int) -> None:
    sam = bam.with_suffix(".sam")
    with open(sam, "wb") as fh:
        cost.run(["minimap2", "-ax", "sr", "-t", str(threads), "-y", ref, str(reads)], stdout=fh)
    with open(bam, "wb") as fh:
        cost.run(["samtools", "sort", "-@", str(threads), "-o", "-", str(sam)], stdout=fh)
    cost.run(["samtools", "index", str(bam)])
    sam.unlink()


def run_umi_tools(sim: dict, out: Path, umi_len: int, threads: int,
                  binary: str) -> tuple[dict[str, str], Cost]:
    cost = Cost()
    tagged = out / "ut_reads.fq.gz"
    _strip_umi_into_name(sim["reads"], tagged, umi_len)
    bam = out / "ut.bam"
    _align(tagged, sim["clones"], bam, cost, threads)
    groups = out / "ut_groups.tsv"
    cost.run([binary, "group", "-I", str(bam), "--group-out", str(groups),
              "--umi-separator", "_", "--log", os.devnull])
    pred: dict[str, str] = {}
    with open(groups) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_read = header.index("read_id")
        i_grp = header.index("final_umi") if "unique_id" not in header else header.index("unique_id")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            pred[f[i_read].rsplit("_", 1)[0]] = f[i_grp]
    return pred, cost


def run_fgbio(sim: dict, out: Path, umi_len: int, threads: int,
              binary: str) -> tuple[dict[str, str], Cost]:
    cost = Cost()
    ubam = out / "fg_unmapped.bam"
    cost.run([binary, "FastqToBam", "-i", sim["reads"], "-o", str(ubam),
              "--read-structures", f"{umi_len}M+T", "--sample", "S1", "--library", "L1"])
    # fgbio wants the UMI carried through the aligner, so the tags go via `samtools fastq -T RX`
    # and minimap2 -y rather than through a re-merge.
    fq = out / "fg_reads.fq"
    with open(fq, "wb") as fh:
        cost.run(["samtools", "fastq", "-T", "RX", str(ubam)], stdout=fh)
    bam = out / "fg.bam"
    _align(fq, sim["clones"], bam, cost, threads)
    grouped = out / "fg_grouped.bam"
    cost.run([binary, "GroupReadsByUmi", "-i", str(bam), "-o", str(grouped),
              "-s", "adjacency", "-e", "1"])
    pred: dict[str, str] = {}
    proc = subprocess.run(["samtools", "view", str(grouped)], capture_output=True, text=True,
                          check=True)
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        mi = next((t.split(":", 2)[2] for t in f[11:] if t.startswith("MI:")), None)
        if mi is not None:
            pred[f[0]] = mi.split("/")[0]
    return pred, cost


TOOLS = ("migec", "umi_tools", "fgbio")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--molecules", type=int, default=20_000)
    ap.add_argument("--clones", type=int, default=200)
    ap.add_argument("--coverage", type=float, default=5.0,
                    help="reads per molecule; 1-3 is the shallow regime everything scales with")
    ap.add_argument("--umi-len", type=int, default=12)
    ap.add_argument("--umi-error", type=float, default=2e-3)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--tools", default=",".join(TOOLS))
    ap.add_argument("--umi-tools-binary", default="umi_tools")
    ap.add_argument("--fgbio-binary", default="fgbio")
    ap.add_argument("--tsv", type=Path, help="also write the table here")
    args = ap.parse_args(argv)

    wanted = [t.strip() for t in args.tools.split(",") if t.strip()]
    for t in wanted:
        if t not in TOOLS:
            raise SystemExit(f"unknown tool {t!r}; pick from {', '.join(TOOLS)}")

    args.out.mkdir(parents=True, exist_ok=True)
    sim = simulate(args.out, args.molecules, args.clones, args.coverage, args.umi_len,
                   args.umi_error, args.seed)
    truth = read_truth(Path(sim["truth_reads"]))
    print(f"# {sim['n_reads']:,} reads, {sim['n_molecules']:,} molecules, "
          f"{sim['n_umi_collisions']:,} UMI collisions, {args.umi_len} nt barcode",
          file=sys.stderr)

    rows: list[dict] = []
    for tool in wanted:
        if tool == "migec":
            pred, cost = run_migec(sim, args.out, args.threads)
        elif tool == "umi_tools":
            if not shutil.which(args.umi_tools_binary):
                print(f"# skipping umi_tools: {args.umi_tools_binary} not on PATH", file=sys.stderr)
                continue
            pred, cost = run_umi_tools(sim, args.out, args.umi_len, max(1, args.threads),
                                       args.umi_tools_binary)
        else:
            if not shutil.which(args.fgbio_binary):
                print(f"# skipping fgbio: {args.fgbio_binary} not on PATH", file=sys.stderr)
                continue
            pred, cost = run_fgbio(sim, args.out, args.umi_len, max(1, args.threads),
                                   args.fgbio_binary)
        scores = adjusted_rand(truth, pred)
        rows.append({"tool": tool, **scores, "seconds": cost.seconds, "peak_rss_bytes": cost.rss})

    cols = ["tool", "reads", "true_molecules", "predicted_clusters", "ari",
            "reads_in_split_molecules", "reads_in_merged_clusters", "seconds", "peak_rss_bytes"]
    lines = ["\t".join(cols)]
    for r in rows:
        cells = [r["tool"]]
        for c in cols[1:]:
            v = r[c]
            cells.append(f"{v:.0f}" if c in ("reads", "true_molecules", "predicted_clusters",
                                             "peak_rss_bytes") else f"{v:.4f}")
        lines.append("\t".join(cells))
    out = "\n".join(lines)
    print(out)
    if args.tsv:
        args.tsv.write_text(out + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
