set terminal svg size 760,520 font "Helvetica,13" enhanced
set output "/Users/mikesh/vcs/code/migec/assets/mig_size_zipf.svg"
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

set title "Molecule size against rank -- a straight line here is Zipf"
set xlabel "rank (molecules at least this deep)"
set ylabel "reads per molecule"
set logscale xy
set key inside bottom left
# The rank curve is the cumulative count of the size spectrum, read from the deep end down, which
# is why the spectrum is emitted at EXACT sizes and not in power-of-two bins: four bins make four
# steps and a straight line cannot be told from a bent one.
zipf = 0
plot "< sort -t'	' -k1,1nr '/Users/mikesh/vcs/code/migec/assets/PBMC.sizes.tsv' | grep -v size"      using (zipf = zipf + $3, zipf):1 with lines lw 2.5 lc rgb "#7570b3"      title "molecules of at least this depth"

