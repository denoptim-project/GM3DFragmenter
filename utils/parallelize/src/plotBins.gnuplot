set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
set xlabel 'MW-bin index'
set ylabel 'number of fragments'
set xtics out rotate 45
set yrange [0:10] 
set grid

plot './all_BINs_Counting_ordered.dat' u 7:xtic(3) ti 'total', '' u 12 ti 'unique '

