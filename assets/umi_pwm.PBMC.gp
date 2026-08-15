set terminal svg size 760,520 font "Helvetica,13" enhanced
set output "/Users/mikesh/vcs/code/migec/assets/umi_pwm.PBMC.svg"
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

set title "Barcode base composition -- PBMC"
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
plot for [i=4:7] "/Users/mikesh/vcs/code/migec/assets/checkout.umi_composition.tsv" using (strcol(1) eq "PBMC" ? $2 : 1/0):i with linespoints      lw 2 pt 7 ps 0.5 lc rgb word("#1b9e77 #d95f02 #7570b3 #e7298a", i - 3)      title word("A C G T", i - 3),      "" using (strcol(1) eq "PBMC" && strcol(3) eq "cell" ? $2 : 1/0):(0)      with points pt 7 ps 0.9 lc rgb "#666666" title "cell barcode positions"

