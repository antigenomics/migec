set terminal svg size 760,520 font "Helvetica,13" enhanced
set output "/Users/mikesh/vcs/code/migec/assets/umi_information.PBMC.svg"
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

set title "What each barcode position is worth -- PBMC"
set xlabel "barcode position (cell barcode, then UMI)"
set ylabel "bits"
set yrange [0:2]
set boxwidth 0.7
set arrow from graph 0, first 2 to graph 1, first 2 nohead dt 2 lc rgb "#999999"
set label "2 bits = a position the synthesiser mixed evenly" at graph 0.02, first 1.9      tc rgb "#999999"
plot "/Users/mikesh/vcs/code/migec/assets/checkout.umi_composition.tsv" using (strcol(1) eq "PBMC" ? $2 : 1/0):9 with boxes lc rgb "#1b9e77"      title "information",      "" using (strcol(1) eq "PBMC" ? $2 : 1/0):8 with linespoints lw 2 pt 7 ps 0.5      lc rgb "#7570b3" title "Shannon entropy"

