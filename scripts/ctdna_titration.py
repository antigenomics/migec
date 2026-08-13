#!/usr/bin/env python3
# 2026-08-13
# What a variant caller is actually given, on a real-UMI ctDNA titration.
#
# PRJNA788522 (Osterlund et al., Clin Chem 2022) is SiMSen-Seq on commercial cell-free DNA
# reference material at KNOWN mutant allele frequencies -- 0% (the `WT` arm), 0.125%, 0.25% and
# 1% -- crossed with DNA input (5/20/80 ng) and sequencing depth (3.3x/10x/30x reads per UMI),
# three replicates each. The UMI is a real 12 nt inline barcode and it survived deposition, which
# is what makes this different from every ctDNA caller benchmark in the literature: those run on
# UMIs simulated after the fact, because the runs they used did not keep theirs.
#
# This script does NOT call variants. It measures the thing that decides whether calling one is
# possible at all: how many MOLECULES survive to the consensus, against how many carry the
# variant at the stated frequency. A caller cannot recover a molecule that consensus never
# produced, so `molecules * VAF` is a hard ceiling on the evidence any caller can see -- and on
# the `WT` arm, where the truth is zero, whatever the pipeline still reports is the floor.
#
# Usage:
#   python scripts/sra_fetch.py get <runs> -o simsen/
#   python scripts/ctdna_titration.py --reads simsen/ --out ctdna/

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

# 20 ng of cell-free DNA is not 20 ng of information. One haploid human genome is 3.3 pg, so a
# 20 ng library holds ~6,060 genome equivalents BEFORE any loss -- and that number, not the read
# depth, is what caps a low-frequency call. Everything below compares against it.
PG_PER_HAPLOID_GENOME = 3.3

# The two studies spell their design into different fields, so there are two grammars and a
# fallback. Never invent the missing pieces: a run whose label says nothing about frequency gets
# `vaf=None` and is reported without a detectability column, rather than being assigned a guess.

# PRJNA788522 sample_alias: <input>ng_<vaf>_<depth>x_rep_<n>. `WT` is the 0% negative control and
# `cell_line` is the undiluted positive, whose frequency the paper does not state.
TITRATION = re.compile(
    r"^(?P<ng>\d+)ng_(?P<vaf>WT|cell_line|[\d.]+)_(?P<depth>[\d.]+)x_rep_(?P<rep>\d+)$")

# PRJNA507366 library_name|experiment_title: `Phusion1|Wildtype`, or
# `3plx Platinum superfi 70ng 0.0625% (3)|0.0625% VAF`.
POLYMERASE = re.compile(r"^(?P<enzyme>[A-Za-z_]+?)(?P<rep>\d*)\|(?P<arm>.*)$")
DILUTION = re.compile(
    r"^.*?(?P<enzyme>Platinum superfi|Plat\w+)\s+(?P<ng>\d+)\s*ng\s+(?P<pct>[\d.]+)%"
    r"(?:\s*\((?P<rep>\d+)\))?\|", re.IGNORECASE)

# `Seracare WT (1) 80ng PlatSuperfi:PlatRegular|Wildtype` -- a two-enzyme comparison at fixed
# input. The enzyme that matters is the one AFTER the colon: the prefix names the library prep,
# the suffix names the polymerase this run actually used.
PAIRED_ENZYME = re.compile(
    r"^Seracare\s+WT\s*\((?P<rep>\d+)\)\s*(?P<ng>\d+)\s*ng\s+\w+:(?P<enzyme>\w+)\|", re.IGNORECASE)


def parse_label(label: str) -> dict:
    """Design fields for one run. Always returns a dict; unknown fields stay empty."""
    base = {"input_ng": "", "arm": label, "vaf": None, "target_depth": "", "replicate": "",
            "enzyme": ""}

    m = TITRATION.match(label)
    if m:
        token = m.group("vaf")
        vaf = None if token == "cell_line" else (0.0 if token == "WT" else float(token) / 100)
        return {**base, "input_ng": int(m.group("ng")), "arm": token, "vaf": vaf,
                "target_depth": float(m.group("depth")), "replicate": int(m.group("rep"))}

    m = DILUTION.match(label)
    if m:
        return {**base, "input_ng": int(m.group("ng")), "arm": f"{m.group('pct')}%",
                "vaf": float(m.group("pct")) / 100, "enzyme": m.group("enzyme"),
                "replicate": int(m.group("rep") or 0)}

    m = PAIRED_ENZYME.match(label)
    if m:
        return {**base, "input_ng": int(m.group("ng")), "arm": "WT", "vaf": 0.0,
                "enzyme": m.group("enzyme"), "replicate": int(m.group("rep"))}

    m = POLYMERASE.match(label)
    if m and m.group("arm").lower().startswith("wildtype"):
        return {**base, "arm": "WT", "vaf": 0.0, "enzyme": m.group("enzyme").rstrip("_"),
                "replicate": int(m.group("rep") or 0)}

    return base


def genome_equivalents(ng: float) -> float:
    return ng * 1000 / PG_PER_HAPLOID_GENOME


def count_amplicons(consensus_fq: Path, prefix_len: int, min_share: float) -> tuple[int, float]:
    """How many amplicons the panel actually has, counted from the consensus rather than assumed.

    Never: a multiplex panel's molecule total is NOT the molecule count at any one site. A variant
    sits on one amplicon, so the number a caller sees is `molecules / amplicons`, and reporting
    the total instead overstates the evidence by exactly the panel size. These are amplicons, so
    every molecule of one starts with the same bases and the count falls out of a prefix tally --
    no reference and no aligner. The long tail of near-miss prefixes is payload error on single
    molecules, which is why there is a share floor.

    Never: the floor cannot be near the noise, because the noise MOVES WITH DEPTH. At 1% this
    counted 5 amplicons on the deepest run and 10 on the shallowest -- not because the panel
    changed but because a shallow consensus carries more payload error, so more error prefixes
    clear the bar. That inflated amplicon count then deflates `molecules_per_amplicon` on exactly
    the runs where the evidence is thinnest, which is backwards. The real panel shows a 7x gap
    between the smallest true amplicon (7.6% of molecules) and the largest error prefix (1.1%),
    so 5% sits in the middle of a gap rather than on a slope, and the count is stable across
    every depth and input in the design.
    """
    import collections
    import gzip

    counts: collections.Counter[str] = collections.Counter()
    total = 0
    with gzip.open(consensus_fq, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                counts[line[:prefix_len]] += 1
                total += 1
    if not total:
        return 0, 0.0
    keep = [v for v in counts.values() if v / total >= min_share]
    return len(keep), sum(keep) / total


def detectable(molecules: float, vaf: float | None, min_support: int) -> dict:
    """The arithmetic a caller inherits.

    `variant_molecules` is the expectation; `p_enough` is the probability that at least
    `min_support` of them are actually there, Poisson in the molecule count. A caller asked to
    find a variant whose expected support is 2 will miss it a third of the time no matter how
    good it is, and no threshold setting recovers a molecule that was never sampled.
    """
    if vaf is None:
        return {"variant_molecules": "", "p_enough": ""}
    lam = molecules * vaf
    # P(X >= min_support) for X ~ Poisson(lam), summed from below because min_support is small.
    below = sum(math.exp(-lam) * lam**k / math.factorial(k) for k in range(min_support))
    return {"variant_molecules": round(lam, 2), "p_enough": round(1 - below, 4)}


def run_one(fq: Path, alias: str, out: Path, umi_len: int, threads: int, rt_floor: float) -> dict:
    """checkout -> refine -> assemble on one run, returning the joined summary.

    `checkout.run` takes a barcode TABLE, not a pattern, and `read_barcodes` does NOT translate a
    slice -- `parse_layout` does, and only the CLI calls it. So a slice written straight into a
    sheet reaches the compiler as the literal string `0:12` and fails on the digit. Translate
    here, the same way the CLI does, rather than hand-writing twelve N's.
    """
    from migec import assemble, checkout, refine
    from migec.cli import _inline_sheet
    from migec.sheet import parse_layout

    sample = fq.name.split(".")[0]
    work = out / sample
    (work / "co").mkdir(parents=True, exist_ok=True)
    pattern, _anchored = parse_layout(f"0:{umi_len}")
    sheet = _inline_sheet(pattern, sample, work / "co")
    co = checkout.run(fq, sheet, work / "co", threads=threads, max_offset=0)
    rf = refine.run(work / "co" / f"{sample}.fq.gz", work / "rf", threads=threads)
    asm = assemble.run(work / "rf" / f"{sample}.fq.gz", work / "as", threads=threads,
                       rt_floor=rt_floor)
    return {"run": sample, "alias": alias, "checkout": co, "refine": rf, "assemble": asm,
            "consensus": work / "as" / f"{sample}.consensus.fq.gz"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reads", type=Path, required=True, help="directory of <run>.fastq.gz")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--design", type=Path,
                    help="TSV of run_accession<TAB>sample_alias; fetched from ENA when absent")
    ap.add_argument("--study", default="PRJNA788522")
    ap.add_argument("--umi-length", type=int, default=12)
    ap.add_argument("--min-support", type=int, default=3,
                    help="molecules a caller needs before it will call (default 3)")
    ap.add_argument("--prefix-length", type=int, default=25,
                    help="consensus prefix used to tell amplicons apart")
    ap.add_argument("--min-share", type=float, default=0.05,
                    help="molecule share an amplicon must hold to be counted; must sit in the gap "
                         "between the smallest amplicon and the largest error prefix, not on a slope")
    ap.add_argument("--rt-floor", type=float, default=1e-4,
                    help="pre-amplification error floor; 1e-4 (Q40) is the one-molecule default")
    ap.add_argument("-t", "--threads", type=int, default=0)
    ap.add_argument("--from-json", action="store_true",
                    help="re-aggregate <out>/raw.json without re-running migec")
    a = ap.parse_args()

    a.out.mkdir(parents=True, exist_ok=True)

    # An optional third column names the study. It matters because the amplicon check is per
    # panel, and two studies in one sheet can legitimately run panels of different sizes.
    design_study: dict[str, str] = {}
    if a.design and a.design.exists():
        design = {}
        for line in a.design.read_text().strip().split("\n")[1:]:
            if not line:
                continue
            parts = line.split("\t")
            design[parts[0]] = parts[1] if len(parts) > 1 else ""
            design_study[parts[0]] = parts[2] if len(parts) > 2 else a.study
    else:
        import urllib.parse
        import urllib.request
        q = urllib.parse.urlencode({"result": "read_run",
                                    "query": f"study_accession={a.study}",
                                    "fields": "run_accession,sample_alias",
                                    "format": "tsv", "limit": 0})
        with urllib.request.urlopen(
                f"https://www.ebi.ac.uk/ena/portal/api/search?{q}", timeout=60) as r:
            rows = r.read().decode().strip().split("\n")[1:]
        design = dict(line.split("\t")[:2] for line in rows if line)
        design_study = dict.fromkeys(design, a.study)

    # The migec pass is ~20 minutes over a hundred runs and the aggregation below is instant, so
    # a change to the summary must not cost a recompute. `raw.json` is written after every run.
    if a.from_json:
        results = json.loads((a.out / "raw.json").read_text())
        for r in results:
            r.update(parse_label(r.get("label", "")))       # re-parse: a grammar may have changed
        print(f"re-aggregating {len(results)} cached runs from {a.out / 'raw.json'}",
              file=sys.stderr)
        return summarise(results, a)

    results = []
    for fq in sorted(a.reads.glob("*.fastq.gz")):
        run = fq.name.split(".")[0].split("_")[0]
        alias = design.get(run, "")
        if not alias:
            print(f"skip {run}: no label in the design table", file=sys.stderr)
            continue
        meta = parse_label(alias)
        print(f"[{len(results) + 1}] {run}  {alias}", file=sys.stderr)
        try:
            r = run_one(fq, alias, a.out, a.umi_length, a.threads, a.rt_floor)
        except Exception as exc:                  # noqa: BLE001 - one bad run must not lose the rest
            print(f"  FAILED {run}: {exc}", file=sys.stderr)
            continue
        rf, asm = r["refine"], r["assemble"]
        reads = rf["reads"]
        mol = asm["molecules"]
        amplicons, amp_share = count_amplicons(r["consensus"], a.prefix_length, a.min_share)
        per_amplicon = mol / amplicons if amplicons else 0
        geq = genome_equivalents(meta["input_ng"]) if meta["input_ng"] != "" else 0
        row = {
            # The raw label is kept so `--from-json` can re-parse it. Storing only the parsed
            # fields would mean a fix to a grammar needs the whole 20-minute pass again.
            "run": run, "study": design_study.get(run, ""), "label": alias, **meta,
            "reads": reads,
            "barcodes": rf["barcodes"],
            "molecules": mol,
            "amplicons": amplicons,
            "amplicon_share": round(amp_share, 4),
            "molecules_per_amplicon": round(per_amplicon),
            "reads_per_molecule": round(reads / mol, 2) if mol else "",
            "barcode_error": round(rf["error_at_depth"], 6),
            "barcode_phred": round(rf["error_phred"], 2),
            "saturated": int(rf["saturated"]),
            "payload_clonality": rf["payload_clonality"],
            "mean_consensus_phred": round(asm["mean_quality"], 2),
            "groups_split": asm["groups_split"],
            "genome_equivalents": round(geq),
            # >1 is expected, not a bug: a double-stranded fragment is barcoded once per STRAND,
            # so the ceiling is 2 per genome equivalent per amplicon before any loss.
            "molecules_per_geq": round(per_amplicon / geq, 2) if geq else "",
            **detectable(per_amplicon, meta["vaf"], a.min_support),
        }
        results.append(row)
        (a.out / "raw.json").write_text(json.dumps(results, indent=1, default=str))

    return summarise(results, a)


def summarise(results: list[dict], a) -> None:  # noqa: ANN001 - argparse.Namespace
    """Everything downstream of the migec pass, so `--from-json` can redo it in a second."""
    if not results:
        raise SystemExit("no runs analysed")

    # The panel size is a property of the STUDY, not of a run. Counting it per run lets the
    # divisor wobble: an undiluted arm carries a mutant amplicon whose prefix splits, and a
    # shallow run loses a weak amplicon below the share floor -- both move the count for reasons
    # that have nothing to do with how many molecules covered a site. Fix it at the study's modal
    # count and keep the per-run count in the table so the disagreement stays visible.
    from statistics import mode
    panel = {s: mode([r["amplicons"] for r in results if r["study"] == s])
             for s in {r["study"] for r in results}}
    for r in results:
        n = panel[r["study"]]
        r["panel_amplicons"] = n
        r["molecules_per_amplicon"] = round(r["molecules"] / n)
        r.update(detectable(r["molecules"] / n, r["vaf"], a.min_support))

    def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
        path.write_text("\t".join(cols) + "\n"
                        + "".join("\t".join(str(r.get(c, "")) for c in cols) + "\n" for r in rows))

    results.sort(key=lambda r: (str(r["enzyme"]), str(r["input_ng"]), str(r["vaf"]),
                                str(r["target_depth"]), str(r["replicate"])))
    write_tsv(a.out / "titration.tsv", list(results[0]), results)

    # A narrower, replicate-averaged table for the docs. Sphinx reads the committed TSV with
    # `csv-table :file:`, so there is no generated .rst to drift from it -- the same rule the
    # figures follow. 100 per-run rows is not a table anyone reads; the per-run detail stays in
    # titration.tsv beside it.
    groups: dict[tuple, list[dict]] = {}
    for r in results:
        groups.setdefault(
            (r["study"], r["enzyme"], r["input_ng"], r["arm"], r["target_depth"]), []).append(r)

    def mean(rows: list[dict], key: str) -> str:
        vals = [r[key] for r in rows if isinstance(r[key], (int, float))]
        return f"{sum(vals) / len(vals):.6g}" if vals else ""

    agg = [{
        "study": k[0], "enzyme": k[1], "input_ng": k[2], "arm": k[3], "depth": k[4],
        "n": len(rows), "amplicons": rows[0]["panel_amplicons"],
        "reads": mean(rows, "reads"),
        "molecules_per_amplicon": mean(rows, "molecules_per_amplicon"),
        "reads_per_molecule": mean(rows, "reads_per_molecule"),
        "barcode_phred": mean(rows, "barcode_phred"),
        "variant_molecules": mean(rows, "variant_molecules"),
        "p_enough": mean(rows, "p_enough"),
    } for k, rows in sorted(groups.items(), key=lambda kv: tuple(str(x) for x in kv[0]))]
    write_tsv(a.out / "titration_summary.tsv", list(agg[0]), agg)
    print(f"\nwrote {a.out / 'titration.tsv'} ({len(results)} runs)", file=sys.stderr)

    # The check on the amplicon threshold: ONE PANEL has one amplicon count, so a count that moves
    # between runs of the same study means --min-share is sitting on the payload-error slope rather
    # than in the gap. Never: this has to be per study -- the two studies here run different
    # panels (PRJNA507366's own labels say "3plx"), so a global check reports a real biological
    # difference as a threshold fault and cries wolf.
    by_study: dict[str, set[int]] = {}
    for r in results:
        by_study.setdefault(r["study"], set()).add(r["amplicons"])
    for study, counts in sorted(by_study.items()):
        n = sum(1 for r in results if r["study"] == study)
        off = sum(1 for r in results
                  if r["study"] == study and r["amplicons"] != panel[study])
        if off:
            print(f"note: {study}: panel is {panel[study]} amplicons; {off}/{n} runs counted "
                  f"something else {sorted(counts)}. The modal size is used for every run -- "
                  f"check those runs if the difference is large.", file=sys.stderr)
        else:
            print(f"{study}: {panel[study]} amplicons, stable across all {n} runs",
                  file=sys.stderr)

    wt = [r for r in results if r["arm"] == "WT"]
    if wt:
        print(f"negative control: {len(wt)} WT runs, "
              f"{sum(r['molecules'] for r in wt):,} molecules with zero true variants",
              file=sys.stderr)
    # A zero-frequency arm has p_enough 0 by definition -- it is the true negative, not a marginal
    # detection. Listing it here buried the rows that actually are marginal under 27 controls.
    marginal = [r for r in results
                if r["p_enough"] != "" and r["vaf"] and r["p_enough"] < 0.9]
    for r in marginal:
        print(f"marginal: {r['run']} {r['arm']} at {r['input_ng']} ng expects "
              f"{r['variant_molecules']} variant molecules, P(>= {a.min_support}) = "
              f"{r['p_enough']}", file=sys.stderr)


if __name__ == "__main__":
    main()
