set terminal svg size 760,520 font "Helvetica,13" enhanced
set output "/Users/mikesh/vcs/code/migec/assets/umi_error_children.svg"
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

set title "Barcode errors against the parent's depth -- distinct children, and the reads in them"
set xlabel "reads carried by the parent"
set ylabel "per parent"
set logscale xy
set key inside top left
# A parent seen c times had c*L barcode bases to miscall, so both series should rise with c -- but
# only one of them can rise forever. There are just 3L distinct barcodes one substitution away, so
# the DISTINCT-children curve bends over and stops; the READS-in-children curve has no ceiling and
# keeps climbing. The dashed line is that ceiling. Where the points leave it is where the barcode
# neighbourhood filled, measured rather than predicted, and it is the same saturation that makes
# the distance-1 error estimate fail downward.
#
# Points, never lines: one row is one exact depth, and past the head most depths hold a handful of
# parents, so a line would draw quantisation noise as structure.
plot "/Users/mikesh/vcs/code/migec/assets/PBMC.umi_errors.tsv" using 1:($5 > 0 ? $5 : 1/0) with points pt 7 ps 0.5 lc rgb "#1b9e77"        title "distinct child barcodes",      "" using 1:($6 > 0 ? $6 : 1/0) with points pt 7 ps 0.5 lc rgb "#d95f02"        title "reads in those children",      "" using 1:7 with lines lw 1.5 dt 2 lc rgb "#808080"        title "3L, the only children there can be"

