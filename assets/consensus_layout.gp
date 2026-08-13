set terminal svg size 760,520 font "Helvetica,13" enhanced
set output "/Users/mikesh/vcs/code/migec/assets/consensus_layout.svg"
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

set title "Contigs per barcode, and how long they came out"
set xlabel "consensus length (nt)"
set ylabel "molecules"
set y2label "barcodes"
set y2tics
set boxwidth 0.9 relative
plot "/Users/mikesh/vcs/code/migec/assets/SRR1763769.mig.tsv" using 8:(1) smooth frequency with boxes lc rgb "#1b9e77" title "length",      "" using 4:(1) axes x1y2 smooth frequency with points pt 7 ps 1.2      lc rgb "#d95f02" title "contigs in the barcode"

