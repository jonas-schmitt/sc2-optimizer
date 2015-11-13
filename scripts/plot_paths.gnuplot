set terminal svg enhanced size 500,500
set output 'paths.svg'
set encoding iso_8859_1
set xrange [0:100]
set yrange [0:100]
unset xtics
unset ytics
set nokey
set style data dots 
set multiplot
call "./scripts/col_counter.gnuplot" './results/pl1_paths_0_0.txt'
plot for [i=1:col_count:2] './results/pl1_paths_0_0.txt' u i:i+1 lc 3 
call "./scripts/col_counter.gnuplot"  './results/pl2_paths_0_0.txt'
plot for [i=1:col_count:2] './results/pl2_paths_0_0.txt' u i:i+1 lc 1 

