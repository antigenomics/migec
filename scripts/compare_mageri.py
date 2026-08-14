#!/usr/bin/env python3
# 2026-08-14
#
# migec 2 against MAGERI 1.1.1, the other descendant of MIGEC 1.
#
# MAGERI is the closer comparator of the two: same author, same UMI model, and it is the tool whose
# assembler this repo's consensus is a rewrite of. It is also a superset -- it assembles, maps and
# calls variants in one run -- so what is compared here is the part they share, the consensus, and
# the extra work MAGERI does is named rather than divided out.
#
# Never: MATCH THE MIG SIZE THRESHOLD. MAGERI's preset carries `forceOverseq=true` and
# `defaultOverseq=5`, so out of the box it assembles only MIGs of 5 reads or more while migec
# defaults to 1. That is the same trap MIGEC 1.2.9's `-m 5` was, and leaving each at its own
# default compares defaults, not algorithms -- and credits migec with recovering molecules MAGERI
# was told to discard. `--min-count` is written into a preset XML and into `--min-reads`, so both
# sides see the same rule.
#
# Note: MAGERI's sub-clustering test is `pcrMinorTestPValue = 0.01`, the nominal threshold X3
# measured as over-calling by 19x on real data (`docs/nulls.rst`). It is left at its default here:
# the point of a head to head is what the other tool does, not what it would do if it were this
# one.
#
# Usage:
#     python scripts/compare_mageri.py --out /tmp/mageri --jar mageri.jar --molecules 20000
#
# Get the jar: gh release download 1.1.1 --repo mikessh/mageri -p mageri.zip  (needs a JDK)

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from compare_migec_v1 import (  # noqa: E402
    ADAPTER,
    _run,
    read_consensus,
    read_truth_consensus,
    score,
)


def write_preset(jar: Path, java: str, out: Path, min_count: int, reads: Path,
                 references: Path, mask: str) -> Path:
    """MAGERI's own preset, with only the MIG size threshold changed.

    Exporting and editing beats writing an XML from scratch: every other parameter stays whatever
    this version of MAGERI ships, which is the thing being compared.
    """
    preset = out / "preset.xml"
    subprocess.run(
        [java, "-jar", str(jar), "--export-preset", str(preset), "-R1", str(reads),
         "--references", str(references), "-M3", mask, "-O", str(out / "_preset_probe")],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
    )
    if not preset.exists():
        raise SystemExit(f"MAGERI wrote no preset to {preset}")
    text = preset.read_text()
    # The count, not the diff: asking for the value it already has is the common case, and a
    # "nothing changed" check would then fail on exactly the run that needed no change.
    edited, n = re.subn(r"<defaultOverseq>\d+</defaultOverseq>",
                        f"<defaultOverseq>{min_count}</defaultOverseq>", text)
    if n != 1:
        raise SystemExit(f"{n} <defaultOverseq> in MAGERI's preset -- the knob has moved")
    preset.write_text(edited)
    return preset


def read_vcf(path: Path, pass_only: bool = True) -> list[dict]:
    """(chrom, pos, ref, alt, af, filter) from a VCF. Both tools put the fraction in INFO/AF."""
    calls = []
    if not path.exists():
        return calls
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 8:
            continue
        flt = f[6]
        if pass_only and flt not in ("PASS", "."):
            continue
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        calls.append({"chrom": f[0], "pos": int(f[1]), "ref": f[3], "alt": f[4],
                      "af": float(info.get("AF", 0.0)), "filter": flt})
    return calls


def read_truth_variants(path: Path) -> dict[tuple[str, int, str], float]:
    truth = {}
    for i, line in enumerate(path.read_text().splitlines()):
        if i == 0 or not line.strip():
            continue
        chrom, pos, _ref, alt, af = line.split("\t")
        truth[(chrom, int(pos), alt)] = float(af)
    return truth


def score_variants(calls: list[dict], truth: dict) -> dict:
    """Sensitivity and false positives against the injected truth.

    A false positive here is a call at a position nothing was injected at, which on a simulated
    library is unambiguous -- the only other source of a non-reference base is sequencing or PCR
    error, and removing that is the entire claim of a UMI pipeline.
    """
    by_key = {(c["chrom"], c["pos"], c["alt"]): c for c in calls}
    hit = set(by_key) & set(truth)
    errs = sorted(abs(by_key[k]["af"] - truth[k]) for k in hit)
    return {
        "truth_variants": len(truth),
        "called": len(calls),
        "tp": len(hit),
        "fp": len(by_key) - len(hit),
        "fn": len(truth) - len(hit),
        "sensitivity": len(hit) / len(truth) if truth else 0.0,
        "precision": len(hit) / len(by_key) if by_key else 0.0,
        "median_af_error": errs[len(errs) // 2] if errs else 0.0,
    }


def mageri_threshold(checkout_txt: Path) -> int:
    """The MIG size threshold MAGERI reports it actually used, not the one we asked for."""
    for line in checkout_txt.read_text().splitlines():
        if line.startswith("#") or line.startswith("sample.group"):
            continue
        return int(line.split("\t")[-1])
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--jar", type=Path, required=True, help="mageri.jar")
    ap.add_argument("--molecules", type=int, default=20_000)
    ap.add_argument("--clones", type=int, default=200)
    ap.add_argument("--coverage", type=float, default=8.0)
    ap.add_argument("--umi-len", type=int, default=12)
    ap.add_argument("--umi-error", type=float, default=2e-3)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--java", default="java")
    ap.add_argument("--min-count", type=int, default=5,
                    help="minimal reads per MIG, applied to BOTH. Defaults to MAGERI's 5 rather "
                         "than migec's 1, because that is the threshold MAGERI was tuned at")
    ap.add_argument("--tsv", type=Path)
    ap.add_argument("--variant-af", type=float, default=0.0,
                    help="turn the clone set into one reference plus point variants at this "
                         "allele fraction, and compare the two pipelines at the VARIANT level. "
                         "0 compares consensuses only")
    ap.add_argument("--caller", default="lofreq", choices=["lofreq", "bcftools"],
                    help="what to run on migec's consensus. MAGERI calls variants itself, so this "
                         "is the half of its pipeline migec does not have")
    ap.add_argument("--tsv-variants", type=Path)
    args = ap.parse_args(argv)

    from tests.synthetic._sim import SimConfig, simulate

    args.out.mkdir(parents=True, exist_ok=True)
    cfg = SimConfig(n_molecules=args.molecules, n_clones=args.clones, coverage=args.coverage,
                    coverage_cv=0.4, umi_len=args.umi_len, umi_error=args.umi_error,
                    adapter=ADAPTER, seed=args.seed, variant_af=args.variant_af)
    sim = simulate(cfg, args.out / "sim")
    truth = read_truth_consensus(Path(sim["truth_consensus"]))
    print(f"# {sim['n_reads']:,} reads, {sim['n_molecules']:,} molecules, "
          f"{len(truth):,} distinct templates, min-count {args.min_count} on both",
          file=sys.stderr)

    rows = []

    # ---------------------------------------------------------------- MAGERI 1.1.1
    # `-M3` is MAGERI's positional rule and it takes the same string migec's pattern grammar does:
    # N for a UMI base, the adapter in lowercase behind it.
    mask = sim["pattern"]
    # With variants there is ONE reference and the variant clones are the answer, so pointing
    # either tool at every clone would let a variant read align to its own reference and be called
    # nothing. Without them the clone set IS the reference set.
    references = sim["reference"] or sim["clones"]
    mag = args.out / "mageri"
    mag.mkdir(parents=True, exist_ok=True)
    preset = write_preset(args.jar, args.java, mag, args.min_count,
                          Path(sim["reads"]), Path(references), mask)
    t, rss = _run([args.java, "-jar", str(args.jar), "-R1", sim["reads"],
                   "--references", references, "-M3", mask, "--import-preset", str(preset),
                   "--sample-name", "S1", "--project-name", "cmp", "-O", str(mag / "out"),
                   "--verbosity", "0"])
    mag_out = mag / "out" / "cmp.S1.assemble.R1.fastq.gz"
    if not mag_out.exists():
        raise SystemExit(f"MAGERI wrote no consensus to {mag_out}")
    used = mageri_threshold(mag / "out" / "cmp.checkout.txt")
    if used != args.min_count:
        raise SystemExit(
            f"MAGERI used a MIG size threshold of {used}, not the {args.min_count} it was given. "
            f"The comparison would be of two different filters"
        )
    rows.append({"tool": "mageri-1.1.1", "seconds": t, "peak_rss_bytes": rss,
                 **score(read_consensus(mag_out), truth)})

    # ---------------------------------------------------------------- migec 2
    from migec.assemble import run as assemble_run
    from migec.checkout import run as checkout_run
    from migec.refine import run as refine_run

    sheet = args.out / "bc.txt"
    sheet.write_text(f"S1\t{mask}\n")
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
    scored = score(read_consensus(v2_out), truth)
    rows.append({"tool": "migec-2", "seconds": seconds,
                 "peak_rss_bytes": st.get("peak_rss_bytes", 0), **scored})

    # MAGERI's clock covers mapping and calling as well, so a third row puts migec on the same
    # footing rather than leaving the reader to discount the first one. The consensus is the same
    # file, so the score is the same; only the clock and the footprint move.
    align_seconds = align_rss = 0.0
    bam = None
    if shutil.which("minimap2") and shutil.which("samtools"):
        bam = v2 / "asm" / "consensus.bam"
        sam = v2 / "asm" / "consensus.sam"
        t3, rss3 = _run(["minimap2", "-ax", "sr", "-y", "-t", str(args.threads or 4),
                         references, str(v2_out), "-o", str(sam)])
        t4, rss4 = _run(["samtools", "sort", "-o", str(bam), str(sam)])
        _run(["samtools", "index", str(bam)])
        align_seconds, align_rss = t3 + t4, max(rss3, rss4)
        rows.append({"tool": "migec-2+minimap2", "seconds": seconds + align_seconds,
                     "peak_rss_bytes": max(st.get("peak_rss_bytes", 0), align_rss), **scored})

    # ---------------------------------------------------------------- the variant level
    # MAGERI calls variants itself; migec stops at the consensus and hands it to a caller. So the
    # comparable thing is the CALL SET, and the caller is named in the row rather than folded in.
    if args.variant_af > 0.0:
        truth_v = read_truth_variants(Path(sim["truth_variants"]))
        vrows = [{"tool": "mageri-1.1.1", "caller": "mageri", "seconds": t,
                  **score_variants(read_vcf(mag / "out" / "cmp.S1.vcf"), truth_v)}]
        if bam is None:
            print("# no minimap2/samtools: migec's call set cannot be produced", file=sys.stderr)
        else:
            vcf = v2 / "asm" / f"{args.caller}.vcf"
            if args.caller == "lofreq":
                tc, _ = _run(["lofreq", "call", "-f", references, "-o", str(vcf), str(bam)])
            else:
                # bcftools is a genotype caller and is here as the fallback, not as a rival: it is
                # documented as calling nothing below ~8% VAF, so a zero from it means the tool.
                pile = v2 / "asm" / "pileup.bcf"
                tc1, _ = _run(["bcftools", "mpileup", "-f", references, "-a", "AD",
                               "-o", str(pile), str(bam)])
                tc2, _ = _run(["bcftools", "call", "-mv", "-o", str(vcf), str(pile)])
                tc = tc1 + tc2
            vrows.append({"tool": "migec-2", "caller": args.caller,
                          "seconds": seconds + align_seconds + tc,
                          **score_variants(read_vcf(vcf), truth_v)})
        # The arm is columns, never a label: one row per (tool, caller, af, coverage) so four runs
        # concatenate into one table that sorts and groups.
        vcols = ["tool", "caller", "variant_af", "coverage", "molecules", "min_count",
                 "truth_variants", "called", "tp", "fp", "fn", "sensitivity", "precision",
                 "median_af_error", "seconds"]
        counts = ("molecules", "min_count", "truth_variants", "called", "tp", "fp", "fn")
        vlines = ["\t".join(vcols)]
        for r in vrows:
            r.update({"variant_af": args.variant_af, "coverage": args.coverage,
                      "molecules": args.molecules, "min_count": args.min_count})
            vlines.append("\t".join(
                str(r[c]) if c in ("tool", "caller") else
                f"{r[c]:.0f}" if c in counts else
                f"{r[c]:g}" if c in ("variant_af", "coverage") else f"{r[c]:.4f}"
                for c in vcols))
        vout = "\n".join(vlines)
        print(vout)
        if args.tsv_variants:
            args.tsv_variants.write_text(vout + "\n")

    cols = ["tool", "min_count", "consensuses", "exact", "precision", "recall", "seconds",
            "peak_rss_bytes"]
    lines = ["\t".join(cols)]
    for r in rows:
        r["min_count"] = args.min_count
        lines.append("\t".join(
            str(r[c]) if c == "tool" else
            f"{r[c]:.0f}" if c in ("min_count", "consensuses", "exact", "peak_rss_bytes")
            else f"{r[c]:.4f}"
            for c in cols))
    out = "\n".join(lines)
    print(out)
    speedup = rows[0]["seconds"] / rows[-1]["seconds"] if rows[-1]["seconds"] else 0.0
    print(f"# migec 2 is {speedup:.1f}x MAGERI's wall clock, against the {rows[-1]['tool']} row. "
          f"Note: MAGERI's clock also covers variant calling, which migec does not do at all",
          file=sys.stderr)
    if args.tsv:
        args.tsv.write_text(out + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
