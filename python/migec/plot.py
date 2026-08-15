"""QC figures, drawn with gnuplot from the tables the pipeline already wrote.

Nothing here computes anything. Every panel is a gnuplot script over a committed TSV, which is the
same rule the stages follow -- a figure must be redrawable from the table next to it, months later,
by someone who does not have the FASTQ any more. That also means `migec plot` needs no plotting
dependency: it writes the ``.gp`` scripts either way, and runs gnuplot when there is one.

Colours are ColorBrewer Dark2 (qualitative, colour-blind safe), never hand-picked, and the four
bases keep the same four colours in every panel that draws them.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

# ColorBrewer Dark2. The first four are A, C, G, T wherever bases are drawn.
DARK2 = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"]

# Mid grey. The one colour that reads on a white page and on a dark one, which is what lets a
# single SVG serve a light README, a dark README and a printed figure without three renders.
INK = "#808080"
# 4:3-ish and not wide. The key lives INSIDE the plot box (see _PREAMBLE), so no horizontal space
# is spent on a legend gutter and the axes get the whole frame.
TERMINALS = {
    # No `background` at all: an SVG with no background rect is transparent, so GitHub's dark mode
    # shows the page behind it rather than a white slab.
    "svg": 'set terminal svg size 760,520 font "Helvetica,13" enhanced',
    "png": 'set terminal pngcairo size 1000,680 font "Helvetica,13" transparent',
    "pdf": 'set terminal pdfcairo size 6.4in,4.4in font "Helvetica,11"',
}

_PREAMBLE = f"""
set output "{{out}}"
set border 3 lw 1 lc rgb "{INK}"
set tics nomirror out textcolor rgb "{INK}"
set title textcolor rgb "{INK}"
set xlabel textcolor rgb "{INK}"
set ylabel textcolor rgb "{INK}"
set y2label textcolor rgb "{INK}"
# Inside, bottom right, opaque-bordered but unfilled -- a legend in the margin makes every figure
# wider than its data and is the first thing a journal asks you to move.
set key inside bottom right box lw 0.5 lc rgb "{INK}" textcolor rgb "{INK}" samplen 2 spacing 1.1
set style fill solid 0.85 noborder
set grid ytics lw 0.5 lc rgb "{INK}" dt 3
set datafile separator "\\t"
set datafile missing "NA"
"""


@dataclass(frozen=True)
class Panel:
    name: str
    source: str  # glob, relative to the directory being plotted
    script: str  # gnuplot; {src} {out} {sample} are substituted
    per_sample: bool = False  # one figure per sample id in column 1


PANELS = [
    # --------------------------------------------------------------------- checkout
    Panel(
        name="umi_pwm",
        source="checkout.umi_composition.tsv",
        per_sample=True,
        script="""
set title "Barcode base composition -- {sample}"
set xlabel "barcode position (cell barcode, then UMI)"
set ylabel "base frequency"
# Never: [0:*], not [0:1]. A well-made barcode sits at 1/4 everywhere, so a fixed unit axis spends
# three quarters of the panel on emptiness and flattens the only thing the figure is for -- the
# departures from 1/4. Anchored at 0 so the bar heights stay honest, and the 1/4 rule stays drawn.
set yrange [0:*]
# One row at the top, over headroom made for it. Five entries stacked at the right sat on the
# curves, because an autoscaled composition fills its own panel and leaves no empty corner.
set key inside top center horizontal maxrows 1
set offsets 0, 0, 0.045, 0
set arrow from graph 0, first 0.25 to graph 1, first 0.25 nohead dt 2 lc rgb "#999999"
set label "1/4" at graph 1.01, first 0.25 tc rgb "#999999"
plot for [i=4:7] "{src}" using (strcol(1) eq "{sample}" ? $2 : 1/0):i with linespoints \
     lw 2 pt 7 ps 0.5 lc rgb word("#1b9e77 #d95f02 #7570b3 #e7298a", i - 3) \
     title word("A C G T", i - 3), \
     "" using (strcol(1) eq "{sample}" && strcol(3) eq "cell" ? $2 : 1/0):(0) \
     with points pt 7 ps 0.9 lc rgb "#666666" title "cell barcode positions"
""",
    ),
    Panel(
        name="umi_information",
        source="checkout.umi_composition.tsv",
        per_sample=True,
        script="""
set title "What each barcode position is worth -- {sample}"
set xlabel "barcode position (cell barcode, then UMI)"
set ylabel "bits"
set yrange [0:2]
set boxwidth 0.7
set arrow from graph 0, first 2 to graph 1, first 2 nohead dt 2 lc rgb "#999999"
set label "2 bits = a position the synthesiser mixed evenly" at graph 0.02, first 1.9 \
     tc rgb "#999999"
plot "{src}" using (strcol(1) eq "{sample}" ? $2 : 1/0):9 with boxes lc rgb "#1b9e77" \
     title "information", \
     "" using (strcol(1) eq "{sample}" ? $2 : 1/0):8 with linespoints lw 2 pt 7 ps 0.5 \
     lc rgb "#7570b3" title "Shannon entropy"
""",
    ),
    Panel(
        name="umi_quality",
        source="checkout.umi_quality.tsv",
        script="""
set title "Reported Phred over the barcode bases"
set xlabel "reported Phred"
set ylabel "bases"
set logscale y
set boxwidth 0.8
plot "{src}" using 2:3 with boxes lc rgb "#1b9e77" title "bases"
""",
    ),
    Panel(
        name="quality_calibration",
        source="checkout.quality_calibration.tsv",
        script="""
set title "What the reported Phred is worth, measured on the pattern's constant bases"
set xlabel "reported Phred"
set ylabel "error probability"
set logscale y
set format y "10^{%T}"
plot "{src}" using 1:5 with lines lw 2 dt 2 lc rgb "#666666" title "nominal 10^{-q/10}", \
     "" using 1:4 with points pt 7 ps 1.2 lc rgb "#d95f02" title "observed", \
     "" using 1:6 with lines lw 2 lc rgb "#1b9e77" title "calibrated"
""",
    ),
    Panel(
        name="coverage",
        source="checkout.coverage.tsv",
        script="""
set title "MIG size distribution -- how many reads each molecule was seen with"
set xlabel "reads per UMI (log2 bins)"
set ylabel "UMIs"
# Base 2 on x because the BINS are powers of two: `checkout.coverage.tsv` is written in MIGEC's
# own doubling bins, so a base-10 axis puts the tick marks somewhere other than the data.
set logscale x 2
set logscale y
set format x "2^{%L}"
set format y "10^{%T}"
# Top left: the distribution falls to the right, so the default bottom-right key sits on the
# deepest bins, which are the ones an over-sequencing question is actually about.
set key inside top left
set boxwidth 0.8 relative
plot "{src}" using 2:4 with boxes lc rgb "#1b9e77" title "UMIs", \
     "" using 2:3 with linespoints lw 2 pt 7 ps 0.6 lc rgb "#d95f02" title "reads in them"
""",
    ),
    Panel(
        name="trimming",
        source="checkout.trimming.tsv",
        script="""
set title "Payload length after trimming -- one spike is a clean trim"
set xlabel "payload length (nt)"
set ylabel "reads"
set logscale y
set boxwidth 0.9 relative
plot "{src}" using 2:3 with boxes lc rgb "#1b9e77" title "reads"
""",
    ),
    Panel(
        name="barcode_space",
        source="checkout.barcode_space.tsv",
        script="""
set title "Barcode space: nominal against what the composition leaves usable"
set ylabel "sequences"
set logscale y
set style data histograms
set style histogram clustered gap 2
set xtics rotate by -30
plot "{src}" using 3:xtic(1) lc rgb "#7570b3" title "nominal 4^L", \
     "" using 4 lc rgb "#1b9e77" title "effective (collision entropy)", \
     "" using 7 lc rgb "#d95f02" title "barcodes observed"
""",
    ),
    # --------------------------------------------------------------------- suggest
    Panel(
        name="cycles",
        source="suggest.cycles.tsv",
        script="""
set title "Per-cycle composition: UMI cycles sit at 1/4, constant cycles at 1"
set xlabel "cycle"
set ylabel "base frequency"
set yrange [0:1]
set y2label "deviation from 1/4"
set y2range [0:0.8]
set y2tics
plot for [i=2:5] "{src}" using 1:i with lines lw 1.5 \
     lc rgb word("#1b9e77 #d95f02 #7570b3 #e7298a", i - 1) title word("A C G T", i - 1), \
     "" using 1:10 axes x1y2 with lines lw 2.5 dt 2 lc rgb "#666666" title "1/4 deviation"
""",
    ),
    Panel(
        name="kmers",
        source="suggest.kmers.tsv",
        script="""
set title "Overrepresented k-mers -- synthetic sequence still in the reads"
set ylabel "observed / expected"
set logscale y
set xtics rotate by -60 font ",10"
set boxwidth 0.7
set key off
plot "< head -21 '{src}'" using 0:4:xtic(1) with boxes lc rgb "#d95f02"
""",
    ),
    # --------------------------------------------------------------------- refine
    Panel(
        name="molecule_rank",
        source="*.rank.tsv",
        script="""
set title "Molecule rank -- reads per molecule, largest first"
set xlabel "molecule rank"
set ylabel "reads"
set logscale xy
set key inside bottom left
plot "{src}" using 1:2 with lines lw 2.5 lc rgb "#1b9e77" title "reads in the molecule"
""",
    ),
    Panel(
        name="cell_rank",
        source="*.cell_rank.tsv",
        script="""
set title "Barcode rank -- unique UMIs per barcode, the knee is where cells stop"
set xlabel "barcode rank"
set ylabel "unique UMIs (molecules)"
set logscale xy
set key inside bottom left
# Cell Ranger's plot, and deliberately the same axes, because it is the figure every user of a
# droplet protocol already knows how to read. UMIs, never reads: one over-amplified molecule would
# otherwise put an empty droplet high on the curve, which is the artefact the plot exists to show.
# The two colours are the call, drawn on the curve rather than described in the caption.
plot "{src}" using 1:($3 == 1 ? $2 : 1/0) with lines lw 3 lc rgb "#1b9e77" title "called cells", \
     "" using 1:($3 == 0 ? $2 : 1/0) with lines lw 2 lc rgb "#999999" title "background", \
     "" using 1:2 with lines lw 0.8 dt 3 lc rgb "{INK}" notitle
""",
    ),
    Panel(
        name="mig_size_spectrum",
        source="*.sizes.tsv",
        script="""
set title "MIG size spectrum -- molecules and the reads they account for"
set xlabel "log(1 + reads per molecule)"
set ylabel "molecules"
set y2label "reads"
set y2tics textcolor rgb "{INK}"
set logscale y
set logscale y2
set boxwidth 0.9 relative
set key inside top right
# Both series, on their own axes, because they peak in different places the moment a library is
# over-sequenced: most MOLECULES are shallow, and most READS are in the deep ones. A figure with
# only the first says the library is fine and a figure with only the second says it is saturated.
# log1p on x so a molecule seen once has a place on the axis; a plain log drops it.
#
# Never: "reads in them" is POINTS, never a line. The spectrum is one row per EXACT size, and past
# the head almost every size holds one molecule -- so reads == size, and a line through those
# points draws the y = x diagonal as the most prominent feature of the figure. It is a tautology,
# not a second mode. Where two or three molecules share a size the same line sawtooths between
# size*1 and size*2, which is integer quantisation drawn as signal, and it bridges gaps in the
# support where no size was observed at all. Points say what the data are: isolated observations,
# one per size, most of them a single molecule.
plot "{src}" using 2:3 with boxes lc rgb "#1b9e77" title "molecules", \
     "" using 2:4 axes x1y2 with points pt 7 ps 0.5 lc rgb "#d95f02" title "reads in them"
""",
    ),
    Panel(
        name="umi_error_children",
        source="*.umi_errors.tsv",
        script="""
set title "Barcode errors against the parent's depth -- distinct children, and the reads in them"
set xlabel "reads carried by the parent"
set ylabel "per parent"
set logscale xy
set key inside top left
# A parent seen c times had c*L barcode bases to miscall, so both series should rise with c -- but
# only one of them can rise forever. There are just 3L distinct barcodes one substitution away, so
# the DISTINCT-children curve bends over and stops; the READS-in-children curve has no ceiling and
# keeps climbing. The dashed line is that ceiling. Where the points leave it is where the barcode
# neighbourhood filled, measured rather than predicted, and it is the same saturation that makes
# the distance-1 error estimate fail downward.
#
# Points, never lines: one row is one exact depth, and past the head most depths hold a handful of
# parents, so a line would draw quantisation noise as structure.
plot "{src}" using 1:($5 > 0 ? $5 : 1/0) with points pt 7 ps 0.5 lc rgb "#1b9e77" \
       title "distinct child barcodes", \
     "" using 1:($6 > 0 ? $6 : 1/0) with points pt 7 ps 0.5 lc rgb "#d95f02" \
       title "reads in those children", \
     "" using 1:7 with lines lw 1.5 dt 2 lc rgb "{INK}" \
       title "3L, the only children there can be"
""",
    ),
    Panel(
        name="umi_error_rate",
        source="*.umi_errors.tsv",
        script="""
set title "Barcode error rate by depth -- two estimators of one number"
set xlabel "reads carried by the parent"
# Phred in the label rather than on a second axis: Phred is -10 log10 of an error rate, so on a log
# y axis one decade IS ten Phred and a linked axis would only relabel the same gridlines. (gnuplot
# refuses to link a nonlinear axis anyway.) It is named because the number this figure produces has
# to be comparable with the barcode's own reported quality, which is the only independent check on
# it there is -- `phred_from_reads` is column 10 of the table if you want it per row.
set ylabel "error per base per read -- one decade is 10 Phred, 1e-3 is Q30"
set logscale xy
set key inside bottom left
# The estimate that matters, and the reason there are two of it. Both invert a model of the same
# eps against the same row:
#
#   distinct children  u(c) = 3L (1 - exp(-c eps / 3))    saturates at 3L
#   reads in children  r(c) = c L eps                     no ceiling
#
# They agree while the neighbourhood is empty and part company as it fills, and the depth where
# they part is worth more than either curve: it is where a distance-1 estimate stops being usable
# on this library. `error_from_variants` is blank past saturation rather than small, because
# inverting a full neighbourhood reports "no errors" for the most error-ridden case there is.
#
# Never: read this at DEPTH. A child whose parent was never sequenced cannot be merged and so
# cannot be counted, which at 1-3 reads/UMI is 80% of them -- the left-hand end of this figure is a
# lower bound and the right-hand end is the measurement.
#
# Never: both series are bounded by the merges correction made, so neither survives a FULL barcode
# space -- there `correct_umis` refuses to merge, rightly, and both fall to zero. Against an
# injected rate: 0.99 and 0.97 of truth at 0.2% occupancy, 0.62 and 0.45 at 33%, and nothing at
# 100%. The `saturated` flag in the report is what says the answer is a floor.
plot "{src}" using 1:($9 > 0 ? $9 : 1/0) with points pt 7 ps 0.5 lc rgb "#d95f02" \
       title "from the reads in the children", \
     "" using 1:($8 > 0 ? $8 : 1/0) with points pt 6 ps 0.5 lc rgb "#1b9e77" \
       title "from the distinct children (saturates)", \
     "" using 1:($11 > 0 ? $11 : 1/0) with lines lw 1.5 dt 2 lc rgb "{INK}" \
       title "what refine reports, read at depth"
""",
    ),
    Panel(
        name="mig_size_zipf",
        source="*.sizes.tsv",
        script="""
set title "Molecule size against rank -- a straight line here is Zipf"
set xlabel "rank (molecules at least this deep)"
set ylabel "reads per molecule"
set logscale xy
set key inside bottom left
# The rank curve is the cumulative count of the size spectrum, read from the deep end down, which
# is why the spectrum is emitted at EXACT sizes and not in power-of-two bins: four bins make four
# steps and a straight line cannot be told from a bent one.
zipf = 0
plot "< sort -t'\t' -k1,1nr '{src}' | grep -v size" \
     using (zipf = zipf + $3, zipf):1 with lines lw 2.5 lc rgb "#7570b3" \
     title "molecules of at least this depth"
""",
    ),
    Panel(
        name="sample_umis",
        source="checkout.summary.tsv",
        script="""
set title "Unique UMIs and reads per sample barcode"
set ylabel "unique UMIs"
set y2label "reads"
set y2tics textcolor rgb "{INK}"
set logscale y
set logscale y2
set style data histograms
set style histogram clustered gap 2
set xtics rotate by -30
set key inside top right
plot "{src}" using 3:xtic(1) lc rgb "#1b9e77" title "unique UMIs", \
     "" using 2 axes x1y2 lc rgb "#d95f02" title "reads"
""",
    ),
    Panel(
        name="refine_coverage",
        source="refine.coverage.tsv",
        script="""
set title "Molecules per MIG size, after barcode correction"
set xlabel "reads per molecule (bin start)"
set ylabel "molecules"
set logscale xy
set boxwidth 0.8 relative
set key off
plot "{src}" using 2:4 with boxes lc rgb "#7570b3"
""",
    ),
    # --------------------------------------------------------------------- assemble
    Panel(
        name="consensus_quality",
        source="assemble.quality_by_depth.tsv",
        script="""
set title "Consensus quality against depth -- the cap is the RT floor, not the instrument"
set xlabel "reads in the molecule (power-of-two bin)"
set ylabel "emitted Phred"
set logscale x 2
set boxwidth 0.6 relative
set key inside bottom right
# A BOX, never a thinned scatter. Emitted quality is discrete and capped at the floor, so at any
# real depth every molecule sits on one or two integers: a cloud of dots draws that as a flat line
# whether the bin holds ten molecules or ten million, and `every 17` then throws away the tails
# that were the only thing the cloud could have shown. These are exact order statistics over every
# molecule, read off the (depth, quality) count grid.
plot "{src}" using 2:6:5:9:8 with candlesticks whiskerbars 0.5 lw 1.5 lc rgb "#1b9e77" \
     title "quartiles and range", \
     "" using 2:7:7:7:7 with candlesticks lw 3 lc rgb "#d95f02" title "median", \
     "" using 2:10 with lines lw 1.5 dt 2 lc rgb "#7570b3" title "mean"
""",
    ),
    Panel(
        name="consensus_error",
        source="*.mig.tsv",
        script="""
set title "Posterior error of the consensus, before the floor is added"
set xlabel "reads in the molecule"
set ylabel "mean posterior error"
set logscale xy
set format y "10^{%T}"
set key inside top right
plot "{src}" using 6:($10 > 0 ? $10 : 1/0) smooth unique with linespoints lw 2.5 pt 7 ps 0.7 \
     lc rgb "#d95f02" title "mean per depth"
""",
    ),
    Panel(
        name="consensus_layout",
        source="*.mig.tsv",
        script="""
set title "Contigs per barcode, and how long they came out"
set xlabel "consensus length (nt)"
set ylabel "molecules"
set y2label "barcodes"
set y2tics
set boxwidth 0.9 relative
plot "{src}" using 8:(1) smooth frequency with boxes lc rgb "#1b9e77" title "length", \
     "" using 4:(1) axes x1y2 smooth frequency with points pt 7 ps 1.2 \
     lc rgb "#d95f02" title "contigs in the barcode"
""",
    ),
    # --------------------------------------------------------------------- benchmarks
    Panel(
        name="benchmark_threads",
        source="benchmark_threads.tsv",
        script="""
set title "checkout thread scaling -- the output is byte-identical at every -t"
set xlabel "threads"
set ylabel "reads / second"
set y2label "peak RSS (MB)"
set y2tics
set logscale x 2
set key inside bottom right
plot "{src}" using 1:3 with linespoints lw 2.5 pt 7 ps 0.8 lc rgb "#1b9e77" \
     title "matching (threads)", \
     "" using 1:2 with linespoints lw 2.5 pt 7 ps 0.8 lc rgb "#d95f02" \
     title "end to end (serial tail included)", \
     "" using 1:4 axes x1y2 with linespoints lw 1.5 dt 2 pt 6 ps 0.8 lc rgb "#7570b3" \
     title "peak RSS"
""",
    ),
    Panel(
        name="assemble_coverage",
        source="assemble.coverage.tsv",
        script="""
set title "Groups per MIG size, as assembled"
set xlabel "reads per group (bin start)"
set ylabel "groups"
set logscale xy
set boxwidth 0.8 relative
set key off
plot "{src}" using 2:4 with boxes lc rgb "#1b9e77"
""",
    ),
]


def _has_rows(path: Path) -> bool:
    """True if the table has a data row under its header. Blank lines do not count."""
    with open(path) as fh:
        next(fh, None)
        return any(line.strip() for line in fh)


def _samples(path: Path) -> list[str]:
    """The distinct values of column 1, in file order. One figure per sample."""
    seen: list[str] = []
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            sample = line.split("\t", 1)[0]
            if sample and sample not in seen:
                seen.append(sample)
    return seen


def run(
    directory: str | Path,
    out_dir: str | Path | None = None,
    fmt: str = "svg",
    gnuplot: str | None = None,
) -> dict:
    """Draw every panel whose table is present in `directory`.

    Returns what was drawn, what was skipped for a missing table, and the gnuplot that ran -- or
    did not. A missing gnuplot is not an error: the scripts are the deliverable and they are
    written either way.
    """
    # Absolute throughout: gnuplot runs with the output directory as its working directory, so a
    # relative path to the table would resolve against the wrong one.
    src_dir = Path(directory).resolve()
    if not src_dir.is_dir():
        raise ValueError(f"{src_dir} is not a directory -- point migec plot at a stage's --out")
    if fmt not in TERMINALS:
        raise ValueError(f"unknown format {fmt!r}; use one of {', '.join(TERMINALS)}")
    out = (Path(out_dir) if out_dir else src_dir / "plots").resolve()
    out.mkdir(parents=True, exist_ok=True)
    # None: find it. A path: use it. "": write the scripts and render nothing, which is what a
    # machine without gnuplot does and therefore what the tests exercise.
    exe = shutil.which("gnuplot") if gnuplot is None else (gnuplot or None)

    drawn: list[str] = []
    skipped: list[str] = []
    empty: list[str] = []
    failed: list[str] = []
    for panel in PANELS:
        matches = sorted(src_dir.glob(panel.source))
        # A table with a header and no rows is the same case as no table at all, and it is a real
        # one rather than a corner: a purely positional pattern (10x) has no constant bases, so
        # `checkout.quality_calibration.tsv` is written empty and there is nothing to calibrate
        # against. Feeding that to gnuplot got "x range is invalid" on stderr and a failure row in
        # the report, which reads as a broken run rather than a chemistry without an anchor.
        present = [m for m in matches if _has_rows(m)]
        if not present:
            # Reported apart, because the advice differs: a missing table means run the stage, an
            # empty one means the stage ran and had nothing to put in it.
            (empty if matches else skipped).append(panel.name)
            continue
        matches = present
        for src in matches:
            stem = panel.name if len(matches) == 1 else f"{panel.name}.{src.name.split('.')[0]}"
            for sample in _samples(src) if panel.per_sample else [""]:
                name = f"{stem}.{sample}" if sample else stem
                figure = out / f"{name}.{fmt}"
                script = (
                    TERMINALS[fmt]
                    + _PREAMBLE.format(out=figure)
                    + panel.script.replace("{src}", str(src))
                    .replace("{sample}", sample)
                    .replace("{INK}", INK)
                )
                gp = out / f"{name}.gp"
                gp.write_text(script.lstrip() + "\n")
                if exe is None:
                    continue
                proc = subprocess.run(
                    [exe, str(gp)], capture_output=True, text=True, cwd=str(out)
                )
                # gnuplot warns on an empty plot rather than failing, so the file is the test.
                if proc.returncode != 0 or not figure.exists():
                    failed.append(f"{name}: {proc.stderr.strip().splitlines()[-1:] or ['no output']}")
                    continue
                drawn.append(figure.name)

    return {
        "out_dir": str(out),
        "format": fmt,
        "gnuplot": exe or "",
        "drawn": drawn,
        "skipped": skipped,
        "empty": empty,
        "failed": failed,
        "scripts": sorted(p.name for p in out.glob("*.gp")),
    }


def format_report(summary: dict) -> str:
    lines = [f"{len(summary['scripts'])} gnuplot scripts in {summary['out_dir']}"]
    if not summary["gnuplot"]:
        lines.append(
            "gnuplot was not found, so nothing was rendered. The scripts are complete and stand "
            "alone: install gnuplot and run `gnuplot *.gp`, or draw the same TSVs anywhere else"
        )
    else:
        lines.append(f"{len(summary['drawn'])} figures drawn with {summary['gnuplot']}")
        for name in summary["drawn"]:
            lines.append(f"  {name}")
    if summary["failed"]:
        lines.append("")
        for f in summary["failed"]:
            lines.append(f"warning: {f}")
    if summary["skipped"]:
        lines.append("")
        lines.append(
            f"no table for: {', '.join(summary['skipped'])} -- run the stage that writes it"
        )
    if summary.get("empty"):
        lines.append("")
        lines.append(
            f"nothing to draw for: {', '.join(summary['empty'])} -- the stage ran and wrote an "
            f"empty table. A purely positional pattern has no constant bases to calibrate against"
        )
    return "\n".join(lines)
