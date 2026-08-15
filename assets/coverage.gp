set terminal svg size 760,520 font "Helvetica,13" enhanced
set output "/Users/mikesh/vcs/code/migec/assets/coverage.svg"
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
plot "/Users/mikesh/vcs/code/migec/assets/checkout.coverage.tsv" using 2:4 with boxes lc rgb "#1b9e77" title "UMIs",      "" using 2:3 with linespoints lw 2 pt 7 ps 0.6 lc rgb "#d95f02" title "reads in them"

