set terminal svg size 900,560 font "Helvetica,13" background rgb "white"
set output "/Users/mikesh/vcs/code/migec/assets/consensus_quality.svg"
set border 3 lw 1 lc rgb "#666666"
set tics nomirror out
set key outside right top box lw 0.5 lc rgb "#cccccc"
set style fill solid 0.85 border lc rgb "#ffffff"
set grid ytics lw 0.5 lc rgb "#e5e5e5"
set datafile separator "\t"
set datafile missing "NA"

set title "Consensus quality against depth -- the cap is the RT floor, not the instrument"
set xlabel "reads in the molecule"
set ylabel "mean emitted Phred"
set logscale x
# `every 17` on the points: a molecule table has millions of rows and an SVG with a million
# circles in it is not a figure. The mean per depth is the claim; the thinned cloud is the spread.
plot "/Users/mikesh/vcs/code/migec/assets/SRR1763769.mig.tsv" using 6:9 every 17 with points pt 7 ps 0.3 lc rgb "#bbbbbb" title "molecules",      "" using 6:9 smooth unique with linespoints lw 2.5 pt 7 ps 0.7 lc rgb "#1b9e77"      title "mean per depth"

