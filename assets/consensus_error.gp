set terminal svg size 900,560 font "Helvetica,13" background rgb "white"
set output "/Users/mikesh/vcs/code/migec/assets/consensus_error.svg"
set border 3 lw 1 lc rgb "#666666"
set tics nomirror out
set key outside right top box lw 0.5 lc rgb "#cccccc"
set style fill solid 0.85 border lc rgb "#ffffff"
set grid ytics lw 0.5 lc rgb "#e5e5e5"
set datafile separator "\t"
set datafile missing "NA"

set title "Posterior error of the consensus, before the floor is added"
set xlabel "reads in the molecule"
set ylabel "mean posterior error"
set logscale xy
set format y "10^{{%%T}}"
plot "/Users/mikesh/vcs/code/migec/assets/SRR1763769.mig.tsv" using 6:($10 > 0 ? $10 : 1/0) every 17 with points pt 7 ps 0.3      lc rgb "#bbbbbb" title "molecules",      "" using 6:($10 > 0 ? $10 : 1/0) smooth unique with linespoints lw 2.5 pt 7 ps 0.7      lc rgb "#d95f02" title "mean per depth"

