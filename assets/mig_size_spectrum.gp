set terminal svg size 760,520 font "Helvetica,13" enhanced
set output "/Users/mikesh/vcs/code/migec/assets/mig_size_spectrum.svg"
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

set title "MIG size spectrum -- molecules and the reads they account for"
set xlabel "log(1 + reads per molecule)"
set ylabel "molecules"
set y2label "reads"
set y2tics textcolor rgb "#808080"
set logscale y
set logscale y2
set boxwidth 0.9 relative
set key inside top right
# Both series, on their own axes, because they peak in different places the moment a library is
# over-sequenced: most MOLECULES are shallow, and most READS are in the deep ones. A figure with
# only the first says the library is fine and a figure with only the second says it is saturated.
# log1p on x so a molecule seen once has a place on the axis; a plain log drops it.
plot "/Users/mikesh/vcs/code/migec/assets/PBMC.sizes.tsv" using 2:3 with boxes lc rgb "#1b9e77" title "molecules",      "" using 2:4 axes x1y2 with lines lw 2.5 lc rgb "#d95f02" title "reads in them"

