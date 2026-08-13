set terminal svg size 900,560 font "Helvetica,13" background rgb "white"
set output "/Users/mikesh/vcs/code/migec/assets/assemble_coverage.svg"
set border 3 lw 1 lc rgb "#666666"
set tics nomirror out
set key outside right top box lw 0.5 lc rgb "#cccccc"
set style fill solid 0.85 border lc rgb "#ffffff"
set grid ytics lw 0.5 lc rgb "#e5e5e5"
set datafile separator "\t"
set datafile missing "NA"

set title "Groups per MIG size, as assembled"
set xlabel "reads per group (bin start)"
set ylabel "groups"
set logscale xy
set boxwidth 0.8 relative
set key off
plot "/Users/mikesh/vcs/code/migec/assets/assemble.coverage.tsv" using 2:4 with boxes lc rgb "#1b9e77"

