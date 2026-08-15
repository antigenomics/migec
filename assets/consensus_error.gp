set terminal svg size 760,520 font "Helvetica,13" enhanced
set output "/Users/mikesh/vcs/code/migec/assets/consensus_error.svg"
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

set title "Posterior error of the consensus, before the floor is added"
set xlabel "reads in the molecule"
set ylabel "mean posterior error"
set logscale xy
set format y "10^{%T}"
set key inside top right
plot "/Users/mikesh/vcs/code/migec/assets/SRR1763769.mig.tsv" using 6:($10 > 0 ? $10 : 1/0) smooth unique with linespoints lw 2.5 pt 7 ps 0.7      lc rgb "#d95f02" title "mean per depth"

