set terminal svg size 760,520 font "Helvetica,13" enhanced
set output "/Users/mikesh/vcs/code/migec/assets/umi_error_rate.svg"
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

set title "Barcode error rate by depth -- two estimators of one number"
set xlabel "reads carried by the parent"
# Phred in the label rather than on a second axis: Phred is -10 log10 of an error rate, so on a log
# y axis one decade IS ten Phred and a linked axis would only relabel the same gridlines. (gnuplot
# refuses to link a nonlinear axis anyway.) It is named because the number this figure produces has
# to be comparable with the barcode's own reported quality, which is the only independent check on
# it there is -- `phred_from_reads` is column 10 of the table if you want it per row.
set ylabel "error per base per read -- one decade is 10 Phred, 1e-3 is Q30"
set logscale xy
set key inside bottom left
# The estimate that matters, and the reason there are two of it. Both invert a model of the same
# eps against the same row:
#
#   distinct children  u(c) = 3L (1 - exp(-c eps / 3))    saturates at 3L
#   reads in children  r(c) = c L eps                     no ceiling
#
# They agree while the neighbourhood is empty and part company as it fills, and the depth where
# they part is worth more than either curve: it is where a distance-1 estimate stops being usable
# on this library. `error_from_variants` is blank past saturation rather than small, because
# inverting a full neighbourhood reports "no errors" for the most error-ridden case there is.
#
# Never: read this at DEPTH. A child whose parent was never sequenced cannot be merged and so
# cannot be counted, which at 1-3 reads/UMI is 80% of them -- the left-hand end of this figure is a
# lower bound and the right-hand end is the measurement.
#
# Never: both series are bounded by the merges correction made, so neither survives a FULL barcode
# space -- there `correct_umis` refuses to merge, rightly, and both fall to zero. Against an
# injected rate: 0.99 and 0.97 of truth at 0.2% occupancy, 0.62 and 0.45 at 33%, and nothing at
# 100%. The `saturated` flag in the report is what says the answer is a floor.
plot "/Users/mikesh/vcs/code/migec/assets/PBMC.umi_errors.tsv" using 1:($9 > 0 ? $9 : 1/0) with points pt 7 ps 0.5 lc rgb "#d95f02"        title "from the reads in the children",      "" using 1:($8 > 0 ? $8 : 1/0) with points pt 6 ps 0.5 lc rgb "#1b9e77"        title "from the distinct children (saturates)",      "" using 1:($11 > 0 ? $11 : 1/0) with lines lw 1.5 dt 2 lc rgb "#808080"        title "what refine reports, read at depth"

