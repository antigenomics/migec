"""The checkout stage: demultiplex, extract UMIs, trim, and write the QC tables.

Everything per-read happens in C++ (`_core.run_checkout`); this module reads the sample sheet,
calls it once, and turns the returned summary into tables. Plots are generated from those tables,
never from inside the C++ -- so a figure can always be redrawn from a committed TSV.
"""

from __future__ import annotations

import json
from pathlib import Path

from migec import _core
from migec.sheet import SampleRow, is_positional, read_barcodes


def run(
    reads: str | Path,
    barcodes: str | Path,
    out_dir: str | Path,
    reads2: str | Path | None = None,
    trim: str = "pattern",
    min_umi_quality: int = 0,
    write_unmatched: bool = False,
    threads: int = 0,
    max_offset: int | None = None,
    limit_reads: int = 0,
    umi_budget_bytes: int = 1 << 30,
    mig: bool = False,
) -> dict:
    """Demultiplex `reads` (and `reads2`, if paired) using `barcodes`, writing into `out_dir`.

    `max_offset=None` picks it: a layout with nothing to score is anchored at the first base,
    anything with an adapter to place it gets a free scan. Never guess the other way -- a free
    scan over an unanchored pattern has no evidence to choose an offset with, and `compile()`
    refuses it rather than picking one.

    `umi_budget_bytes` caps what the UMI counters may hold resident for the whole run. Past it
    they range-partition into `<out_dir>/.umi_spill`, which is removed when the summary is
    written, and every statistic over them streams one bucket at a time. It is not a CLI flag: 1
    GB is the point past which a counter is a problem on any machine, and the stage below it costs
    nothing. 0 keeps everything resident, which is what the counters did before they could
    partition.

    `mig=True` writes `<sample>.<bbb>.mig` buckets instead of one FASTQ per sample. They come out
    range-partitioned on the key `assemble` groups by, so `assemble` reads them directly and skips
    its own partition pass. FASTQ is the default and stays it: it is what every aligner and every
    existing pipeline speaks, and a `.mig` file is an intermediate only migec reads.
    """
    rows: list[SampleRow] = read_barcodes(barcodes)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    if max_offset is None:
        specs = [r.pattern for r in rows] + [r.slave for r in rows if r.slave]
        max_offset = 0 if all(is_positional(s) for s in specs) else -1

    summary = _core.run_checkout(
        str(reads),
        "" if reads2 is None else str(reads2),
        [r.sample_id for r in rows],
        [r.pattern for r in rows],
        [r.slave or "" for r in rows],
        str(out) + "/",
        trim,
        min_umi_quality,
        write_unmatched,
        threads,
        max_offset,
        limit_reads,
        umi_budget_bytes,
        mig,
    )
    summary["input"] = str(reads)
    summary["input2"] = "" if reads2 is None else str(reads2)
    summary["barcodes"] = str(barcodes)
    summary["patterns"] = {r.sample_id: r.pattern for r in rows}
    summary["slaves"] = {r.sample_id: r.slave for r in rows if r.slave}

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

    # The birthday arithmetic and the error budget: every number here is derived from the two
    # above, but they are the numbers that say whether a molecule count means anything.
    with open(out / "checkout.barcode_space.tsv", "w") as fh:
        fh.write(
            "sample_id\tumi_length\tnominal_space\teffective_space\teffective_length\tbias_loss\t"
            "observed_barcodes\toccupancy\tlambda\tmolecules\thidden\tp_multi\tsaturated\t"
            "err_from_phred\tmean_phred\terr_from_polymerase\terr_predicted\terr_estimated\t"
            "err_ratio\tbarcodes_with_error\tneighbour_occupancy\terr_unreliable\n"
        )
        for x in summary["samples"]:
            b, e = x["barcode_space"], x["error_budget"]
            fh.write(
                f"{x['sample_id']}\t{b['length']}\t{b['nominal_space']:.0f}\t"
                f"{b['effective_space']:.1f}\t{b['effective_length']:.4f}\t{b['bias_loss']:.6f}\t"
                f"{b['observed']}\t{b['occupancy']:.6f}\t{b['lambda']:.6f}\t"
                f"{b['molecules']:.1f}\t{b['hidden']:.1f}\t{b['p_multi']:.6f}\t"
                f"{int(b['saturated'])}\t{e['from_phred']:.6e}\t{e['mean_phred']:.2f}\t"
                f"{e['from_polymerase']:.6e}\t{e['predicted']:.6e}\t{e['estimated']:.6e}\t"
                f"{e['ratio']:.4f}\t{e['barcodes_with_error']:.6f}\t"
                f"{e['neighbour_occupancy']:.6f}\t{int(e['estimate_unreliable'])}\n"
            )

    # What the reported Phred is worth, measured against the pattern's own constant bases.
    with open(out / "checkout.quality_calibration.tsv", "w") as fh:
        fh.write("phred\tbases\tmismatches\tobserved\tnominal\tcalibrated\n")
        for q in summary["quality_calibration"]["per_phred"]:
            fh.write(
                f"{q['phred']}\t{q['bases']}\t{q['mismatches']}\t{q['observed']:.6e}\t"
                f"{q['nominal']:.6e}\t{q['fitted']:.6e}\n"
            )
    with open(out / "checkout.pattern_positions.tsv", "w") as fh:
        fh.write("position\tbases\tmismatches\trate\tused_for_calibration\n")
        for p in summary["quality_calibration"]["per_position"]:
            fh.write(
                f"{p['position']}\t{p['bases']}\t{p['mismatches']}\t{p['rate']:.6e}\t"
                f"{int(p['used'])}\n"
            )

    # Reported Phred over the barcode bases -- the input to the predicted error rate.
    with open(out / "checkout.umi_quality.tsv", "w") as fh:
        fh.write("sample_id\tphred\tbases\terror_probability\n")
        for x in summary["samples"]:
            for q in x["umi_phred"]:
                fh.write(f"{x['sample_id']}\t{q['phred']}\t{q['bases']}\t{10 ** (-q['phred'] / 10):.6e}\n")

    # What the trim left: the payload length distribution per sample. A pattern matched one base
    # off still matches, and this is where that shows.
    with open(out / "checkout.trimming.tsv", "w") as fh:
        fh.write("sample_id\tpayload_length\treads\tat_least\n")
        for s in summary["samples"]:
            for L in s["payload_lengths"]:
                fh.write(
                    f"{s['sample_id']}\t{L['length']}\t{L['reads']}\t{int(L['at_least'])}\n"
                )

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
    # Two clocks, because they scale differently: matching is what `-t` speeds up, and the UMI
    # statistics are a pass over every distinct barcode that runs once per sample at the end. Their
    # distance-1 census threads now, so the tail is no longer wholly serial -- but it is still a
    # separate wall, and reporting one number would publish a throughput no user sees.
    total, match = c.get("wall_seconds", 0.0), c.get("match_seconds", 0.0)
    lines.append(
        f"{_dur(total)} "
        f"({c.get('reads_per_second', 0.0):,.0f} reads/s) = "
        f"{_dur(match)} matching on {c.get('threads', 1)} threads "
        f"+ {_dur(max(0.0, total - match))} UMI statistics"
    )
    lines.append(
        f"peak RSS {_bytes(c.get('peak_rss_bytes', 0))} "
        f"of which UMI counters {_bytes(c.get('umi_memory_bytes', 0))}"
        + (" (range-partitioned to disk)" if c.get("umi_spilled") else "")
    )
    lines.append("")
    lines.append(
        f"{'sample':<12}{'reads':>12}{'UMIs':>12}{'reads/UMI':>11}{'UMI len':>9}{'eff len':>9}"
        f"{'payload':>9}"
    )
    for s in summary["samples"]:
        lines.append(
            f"{s['sample_id']:<12}{s['reads']:>12,}{s['umis']:>12,}"
            f"{s['mean_reads_per_umi']:>11.2f}{s['umi_length']:>9}{s['effective_length']:>9.2f}"
            f"{s['mean_payload_length']:>9.1f}"
        )

    # The birthday arithmetic, per sample. Occupancy is the number that decides whether the
    # molecule count means anything, and it is not visible from reads/UMI alone.
    lines.append("")
    lines.append(
        f"{'sample':<12}{'space':>13}{'occupancy':>11}{'MIGs >1 mol':>13}{'molecules':>12}"
        f"{'err pred':>11}{'err est':>10}"
    )
    for s in summary["samples"]:
        b, e = s["barcode_space"], s["error_budget"]
        lines.append(
            f"{s['sample_id']:<12}{b['effective_space']:>13,.0f}{100 * b['occupancy']:>10.1f}%"
            f"{100 * b['p_multi']:>12.1f}%{b['molecules']:>12,.0f}"
            f"{e['predicted']:>11.1e}{e['estimated']:>10.1e}"
        )

    cal = c.get("quality_calibration", {})
    if cal.get("fitted"):
        lines.append("")
        lines.append(
            f"reported Phred is worth {cal['slope']:.2f}x its nominal error, measured on "
            f"{cal['bases']:,} constant pattern bases"
        )
        lines.append(
            f"  the fit's intercept is {cal['quality_independent']:.1e} per base -- the "
            f"SYNTHESISED anchor's own defect rate, not a sequencing floor"
        )
        if cal["positions_dropped"]:
            lines.append(
                f"  {cal['positions_dropped']} pattern position(s) dropped as variable rather "
                f"than miscalled"
            )

    warnings = []
    # The UMI counters were the one allocation that grew with the library rather than with the
    # chunk size, so they were the thing that decided whether a run fits. They partition now, and
    # a counter that is still over the budget is one that could not: a barcode too short to hold
    # two partition prefixes, or spilling switched off.
    if c.get("umi_memory_bytes", 0) > 1 << 30 and not c.get("umi_spilled"):
        warnings.append(
            f"UMI counters hold {_bytes(c['umi_memory_bytes'])} and did not partition. This grows "
            f"with the number of distinct UMIs, so a much larger input may not fit in memory -- "
            f"check that umi_budget_bytes is not 0 and that the barcode is long enough to split"
        )
    if cal.get("fitted") and cal["slope"] > 2.0:
        warnings.append(
            f"the reported Phred understates the error by {cal['slope']:.1f}x. Every likelihood "
            f"here is computed from it, so the calibrated table is the one to use downstream"
        )
    if cal.get("fitted") and cal["quality_independent"] > 1e-3:
        warnings.append(
            f"the constant pattern bases mismatch at {cal['quality_independent']:.1e} per base "
            f"even at the best quality. That is the oligo, not the instrument: synthesis runs "
            f"about one defect per 200-500 bases. It caps how well any barcode on this primer "
            f"can be read"
        )
    # Zero is not "a low rate", it is a configuration error, and it is the one outcome a summary
    # table renders indistinguishably from a successful run of an empty file. Say what it means.
    if c["total"] and c["assigned"] == 0:
        warnings.append(
            "NOTHING matched. The pattern was never placed in any read, which is a declaration "
            "error rather than a data one. Check with `migec suggest`, and note that an anchor "
            "shorter than four bases cannot be scored: three matching bases are 6.0 bits against "
            "an anchored acceptance bar of 6.64, so a pattern ending in a short constant run "
            "refuses every read. Write those bases as `.` (skipped) instead"
        )
    elif c["total"] and c["assigned"] / c["total"] < 0.5:
        warnings.append(
            "less than half of reads matched a pattern -- run `migec suggest` to check where the "
            "barcode actually is"
        )
    # Every warning below is computed from the reads that matched. With none, they are all
    # divisions by zero dressed as findings -- "base composition costs 100% of the barcode space"
    # is not a fact about the library, and printing it buries the one line that matters.
    if c["total"] and c["assigned"] == 0:
        lines.append("")
        lines.append(f"warning: {warnings[-1]}")
        return "\n".join(lines)

    # One line however many samples are under-sequenced: on a 96-plex sheet the per-sample form
    # buries every other warning.
    thin = [s for s in summary["samples"] if not s["over_sequenced"]]
    if thin:
        worst = min(thin, key=lambda s: s["mean_reads_per_umi"])
        which = thin[0]["sample_id"] if len(thin) == 1 else f"{len(thin)} samples"
        warnings.append(
            f"{which} under-sequenced (as low as {worst['mean_reads_per_umi']:.1f} reads/UMI in "
            f"{worst['sample_id']}). Consensus assembly needs over-sequencing; below ~5 most "
            f"molecules are seen once"
        )
    for s in summary["samples"]:
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
        b, e = s["barcode_space"], s["error_budget"]
        # The birthday problem is not a warning until it is: a barcode space that is 30% full puts
        # a sixth of all MIGs on two or more molecules, and no consensus repairs that.
        if b["p_multi"] > 0.05:
            warnings.append(
                f"{s['sample_id']}: {100 * b['p_multi']:.0f}% of MIGs hold more than one molecule "
                f"({100 * b['occupancy']:.0f}% of a {b['effective_space']:,.0f} barcode space is "
                f"occupied). Their consensus is a mixture of templates, not a molecule"
            )
        if b["bias_loss"] > 0.25:
            warnings.append(
                f"{s['sample_id']}: base composition costs {100 * b['bias_loss']:.0f}% of the "
                f"barcode space -- 4^{b['length']} = {b['nominal_space']:,.0f} nominal against "
                f"{b['effective_space']:,.0f} usable. That is the synthesiser mix, and it makes "
                f"collisions more frequent than the length suggests"
            )
        if e["estimate_unreliable"]:
            warnings.append(
                f"{s['sample_id']}: the barcode error estimate ({e['estimated']:.1e}) is not "
                f"reliable here -- {100 * e['neighbour_occupancy']:.0f}% of each barcode's "
                f"1-substitution neighbourhood is itself occupied, so the distance-1 excess it "
                f"reads is a small difference of two large numbers. Phred and polymerase predict "
                f"{e['predicted']:.1e}"
            )
        elif e["predicted"] > 0 and not 0.3 < e["ratio"] < 3.0:
            warnings.append(
                f"{s['sample_id']}: barcode error estimated at {e['estimated']:.1e} against "
                f"{e['predicted']:.1e} predicted from the reported Phred ({e['from_phred']:.1e}) "
                f"and polymerase ({e['from_polymerase']:.1e}) -- {e['ratio']:.2f}x. One of the two "
                f"is wrong; on a 2-colour instrument suspect the nominal Phred first"
            )
    if warnings:
        lines.append("")
        for w in warnings:
            lines.append(f"warning: {w}")
    return "\n".join(lines)


def _pct(n: int, total: int) -> str:
    return f"{100.0 * n / total:.1f}%" if total else "n/a"


def _bytes(n: float) -> str:
    for unit in ("B", "kB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    raise AssertionError("unreachable: the loop returns at TB")


def _dur(s: float) -> str:
    if s < 60:
        return f"{s:.1f} s"
    return f"{int(s // 60)}m{s % 60:04.1f}s"
