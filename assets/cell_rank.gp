set terminal svg size 760,520 font "Helvetica,13" enhanced
set output "/Users/mikesh/vcs/code/migec/assets/cell_rank.svg"
set border 3 lw 1 lc rgb "#808080"
set tics nomirror out textcolor rgb "#808080"
set title textcolor rgb "#808080"
set xlabel textcolor rgb "#808080"
set ylabel textcolor rgb "#808080"
set y2label textcolor rgb "#808080"
# Inside, bottom right, opaque-bordered but unfilled -- a legend in the margin makes every figure
# wider than its data and is the first thing a journal asks you to move.
set key inside bottom right box lw 0.5 lc rgb "#808080" textcolor rgb "#808080" samplen 2 spacing 1.1
set style fill solid 0.85 noborder
set grid ytics lw 0.5 lc rgb "#808080" dt 3
set datafile separator "\t"
set datafile missing "NA"

set title "Barcode rank -- unique UMIs per barcode, the knee is where cells stop"
set xlabel "barcode rank"
set ylabel "unique UMIs (molecules)"
set logscale xy
set key inside bottom left
# Cell Ranger's plot, and deliberately the same axes, because it is the figure every user of a
# droplet protocol already knows how to read. UMIs, never reads: one over-amplified molecule would
# otherwise put an empty droplet high on the curve, which is the artefact the plot exists to show.
# The two colours are the call, drawn on the curve rather than described in the caption.
plot "/Users/mikesh/vcs/code/migec/assets/PBMC.cell_rank.tsv" using 1:($3 == 1 ? $2 : 1/0) with lines lw 3 lc rgb "#1b9e77" title "called cells",      "" using 1:($3 == 0 ? $2 : 1/0) with lines lw 2 lc rgb "#999999" title "background",      "" using 1:2 with lines lw 0.8 dt 3 lc rgb "#808080" notitle

