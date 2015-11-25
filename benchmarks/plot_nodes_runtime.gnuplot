set terminal svg enhanced size 500,500
set output 'nodes_runtime.svg'
set encoding iso_8859_1
set xlabel 'Sockets'
set ylabel 'Runtime (s)'
set xrange [4:32]
set yrange [0:1200]
set xtics 4,4,32 
set nokey
set style data line
plot "weak_scaling_nodes.out" u 1:($2/1000.0)
