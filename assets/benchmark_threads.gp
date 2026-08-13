set terminal svg size 900,560 font "Helvetica,13" background rgb "white"
set output "/Users/mikesh/vcs/code/migec/assets/benchmark_threads.svg"
set border 3 lw 1 lc rgb "#666666"
set tics nomirror out
set key outside right top box lw 0.5 lc rgb "#cccccc"
set style fill solid 0.85 border lc rgb "#ffffff"
set grid ytics lw 0.5 lc rgb "#e5e5e5"
set datafile separator "\t"
set datafile missing "NA"

set title "checkout thread scaling -- the output is byte-identical at every -t"
set xlabel "threads"
set ylabel "reads / second"
set y2label "peak RSS (MB)"
set y2tics
set logscale x 2
set key inside bottom right
plot "/Users/mikesh/vcs/code/migec/assets/benchmark_threads.tsv" using 1:3 with linespoints lw 2.5 pt 7 ps 0.8 lc rgb "#1b9e77"      title "matching (threads)",      "" using 1:2 with linespoints lw 2.5 pt 7 ps 0.8 lc rgb "#d95f02"      title "end to end (serial tail included)",      "" using 1:4 axes x1y2 with linespoints lw 1.5 dt 2 pt 6 ps 0.8 lc rgb "#7570b3"      title "peak RSS"

