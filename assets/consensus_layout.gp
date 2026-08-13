set terminal svg size 900,560 font "Helvetica,13" background rgb "white"
set output "/Users/mikesh/vcs/code/migec/assets/consensus_layout.svg"
set border 3 lw 1 lc rgb "#666666"
set tics nomirror out
set key outside right top box lw 0.5 lc rgb "#cccccc"
set style fill solid 0.85 border lc rgb "#ffffff"
set grid ytics lw 0.5 lc rgb "#e5e5e5"
set datafile separator "\t"
set datafile missing "NA"

set title "Contigs per barcode, and how long they came out"
set xlabel "consensus length (nt)"
set ylabel "molecules"
set y2label "barcodes"
set y2tics
set boxwidth 0.9 relative
plot "/Users/mikesh/vcs/code/migec/assets/SRR1763769.mig.tsv" using 8:(1) smooth frequency with boxes lc rgb "#1b9e77" title "length",      "" using 4:(1) axes x1y2 smooth frequency with points pt 7 ps 1.2      lc rgb "#d95f02" title "contigs in the barcode"

