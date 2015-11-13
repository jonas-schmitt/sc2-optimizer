set terminal svg enhanced size 500,500
set output 'paths.svg'
set encoding iso_8859_1
set xrange [0:100]
set yrange [0:100]
set nokey
set style data dots 
plot "pl1_paths_0_0.txt" u 1:2 lc 3, \
"pl1_paths_0_0.txt" u 3:4 lc 3, \
"pl1_paths_0_0.txt" u 5:6 lc 3, \
"pl1_paths_0_0.txt" u 7:8 lc 3, \
"pl1_paths_0_0.txt" u 9:10 lc 3, \
"pl1_paths_0_0.txt" u 11:12 lc 3, \
"pl2_paths_0_0.txt" u 1:2 lc 1, \
"pl2_paths_0_0.txt" u 3:4 lc 1, \
"pl2_paths_0_0.txt" u 5:6 lc 1, \
"pl2_paths_0_0.txt" u 7:8 lc 1
"pl2_paths_0_0.txt" u 9:10 lc 1, \
"pl2_paths_0_0.txt" u 11:12 lc 1
