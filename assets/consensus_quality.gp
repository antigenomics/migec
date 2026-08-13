set terminal svg size 760,520 font "Helvetica,13" enhanced
set output "/Users/mikesh/vcs/code/migec/assets/consensus_quality.svg"
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

set title "Consensus quality against depth -- the cap is the RT floor, not the instrument"
set xlabel "reads in the molecule (power-of-two bin)"
set ylabel "emitted Phred"
set logscale x 2
set boxwidth 0.6 relative
set key inside bottom right
# A BOX, never a thinned scatter. Emitted quality is discrete and capped at the floor, so at any
# real depth every molecule sits on one or two integers: a cloud of dots draws that as a flat line
# whether the bin holds ten molecules or ten million, and `every 17` then throws away the tails
# that were the only thing the cloud could have shown. These are exact order statistics over every
# molecule, read off the (depth, quality) count grid.
plot "/Users/mikesh/vcs/code/migec/assets/assemble.quality_by_depth.tsv" using 2:6:5:9:8 with candlesticks whiskerbars 0.5 lw 1.5 lc rgb "#1b9e77"      title "quartiles and range",      "" using 2:7:7:7:7 with candlesticks lw 3 lc rgb "#d95f02" title "median",      "" using 2:10 with lines lw 1.5 dt 2 lc rgb "#7570b3" title "mean"

