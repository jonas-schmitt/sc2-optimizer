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
#call "./scripts/col_counter.gnuplot" './results/pl1_paths_0_2.txt'
#plot for [i=1:col_count:2] './results/pl1_paths_0_2.txt' u i:i+1 lc 3 
#call "./scripts/col_counter.gnuplot"  './results/pl2_paths_0_2.txt'
#plot for [i=1:col_count:2] './results/pl2_paths_0_2.txt' u i:i+1 lc 1 

plot "./results/pl1_paths_0_2.txt" u 1:2 lc rgb "dark-red", \
"./results/pl1_paths_0_2.txt" u 3:4 lc rgb "dark-red", \
"./results/pl1_paths_0_2.txt" u 5:6 lc rgb "dark-red", \
"./results/pl1_paths_0_2.txt" u 7:8 lc rgb "dark-red", \
"./results/pl1_paths_0_2.txt" u 9:10 lc rgb "dark-red", \
"./results/pl1_paths_0_2.txt" u 11:12 lc rgb "red", \
"./results/pl1_paths_0_2.txt" u 13:14 lc rgb "red", \
"./results/pl2_paths_0_2.txt" u 1:2 lc rgb "blue", \
"./results/pl2_paths_0_2.txt" u 3:4 lc rgb "blue", \
"./results/pl2_paths_0_2.txt" u 5:6 lc rgb "blue", \
"./results/pl2_paths_0_2.txt" u 7:8 lc rgb "blue", \
"./results/pl2_paths_0_2.txt" u 9:10 lc rgb "blue", \
"./results/pl2_paths_0_2.txt" u 11:12 lc rgb "dark-blue", \
"./results/pl2_paths_0_2.txt" u 13:14 lc rgb "dark-blue"
