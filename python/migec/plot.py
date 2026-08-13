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

TERMINALS = {
    "svg": 'set terminal svg size 900,560 font "Helvetica,13" background rgb "white"',
    "png": 'set terminal pngcairo size 1100,680 font "Helvetica,13"',
    "pdf": 'set terminal pdfcairo size 7in,4.4in font "Helvetica,11"',
}

_PREAMBLE = """
set output "{out}"
set border 3 lw 1 lc rgb "#666666"
set tics nomirror out
set key outside right top box lw 0.5 lc rgb "#cccccc"
set style fill solid 0.85 border lc rgb "#ffffff"
set grid ytics lw 0.5 lc rgb "#e5e5e5"
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
set title "UMI base composition -- {sample}"
set xlabel "barcode position"
set ylabel "base frequency"
set yrange [0:1]
set arrow from graph 0, first 0.25 to graph 1, first 0.25 nohead dt 2 lc rgb "#999999"
set label "1/4" at graph 1.01, first 0.25 tc rgb "#999999"
plot for [i=3:6] "{src}" using (strcol(1) eq "{sample}" ? $2 : 1/0):i with linespoints \
     lw 2 pt 7 ps 0.5 lc rgb word("#1b9e77 #d95f02 #7570b3 #e7298a", i - 2) \
     title word("A C G T", i - 2)
""",
    ),
    Panel(
        name="umi_information",
        source="checkout.umi_composition.tsv",
        per_sample=True,
        script="""
set title "What each barcode position is worth -- {sample}"
set xlabel "barcode position"
set ylabel "bits"
set yrange [0:2]
set boxwidth 0.7
set arrow from graph 0, first 2 to graph 1, first 2 nohead dt 2 lc rgb "#999999"
set label "2 bits = a position the synthesiser mixed evenly" at graph 0.02, first 1.9 \
     tc rgb "#999999"
plot "{src}" using (strcol(1) eq "{sample}" ? $2 : 1/0):8 with boxes lc rgb "#1b9e77" \
     title "information", \
     "" using (strcol(1) eq "{sample}" ? $2 : 1/0):7 with linespoints lw 2 pt 7 ps 0.5 \
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
set format y "10^{{%%T}}"
plot "{src}" using 1:5 with lines lw 2 dt 2 lc rgb "#666666" title "nominal 10^{{-q/10}}", \
     "" using 1:4 with points pt 7 ps 1.2 lc rgb "#d95f02" title "observed", \
     "" using 1:6 with lines lw 2 lc rgb "#1b9e77" title "calibrated"
""",
    ),
    Panel(
        name="coverage",
        source="checkout.coverage.tsv",
        script="""
set title "MIG size distribution -- how many reads each molecule was seen with"
set xlabel "reads per UMI"
set ylabel "UMIs"
set logscale xy
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
        name="cell_rank",
        source="*.rank.tsv",
        script="""
set title "Barcode rank -- the knee is where cells stop and ambient starts"
set xlabel "barcode rank"
set ylabel "reads"
set logscale xy
set key off
plot "{src}" using 1:2 with lines lw 2.5 lc rgb "#1b9e77"
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
        source="*.mig.tsv",
        script="""
set title "Consensus quality against depth -- the cap is the RT floor, not the instrument"
set xlabel "reads in the molecule"
set ylabel "mean emitted Phred"
set logscale x
# `every 17` on the points: a molecule table has millions of rows and an SVG with a million
# circles in it is not a figure. The mean per depth is the claim; the thinned cloud is the spread.
plot "{src}" using 6:9 every 17 with points pt 7 ps 0.3 lc rgb "#bbbbbb" title "molecules", \
     "" using 6:9 smooth unique with linespoints lw 2.5 pt 7 ps 0.7 lc rgb "#1b9e77" \
     title "mean per depth"
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
set format y "10^{{%%T}}"
plot "{src}" using 6:($10 > 0 ? $10 : 1/0) every 17 with points pt 7 ps 0.3 \
     lc rgb "#bbbbbb" title "molecules", \
     "" using 6:($10 > 0 ? $10 : 1/0) smooth unique with linespoints lw 2.5 pt 7 ps 0.7 \
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
    failed: list[str] = []
    for panel in PANELS:
        matches = sorted(src_dir.glob(panel.source))
        if not matches:
            skipped.append(panel.name)
            continue
        for src in matches:
            stem = panel.name if len(matches) == 1 else f"{panel.name}.{src.name.split('.')[0]}"
            for sample in _samples(src) if panel.per_sample else [""]:
                name = f"{stem}.{sample}" if sample else stem
                figure = out / f"{name}.{fmt}"
                script = (
                    TERMINALS[fmt]
                    + _PREAMBLE.format(out=figure)
                    + panel.script.replace("{src}", str(src)).replace("{sample}", sample)
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
    return "\n".join(lines)
