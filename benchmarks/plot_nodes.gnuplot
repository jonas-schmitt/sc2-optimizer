set terminal svg enhanced size 500,500
set output 'weak_scaling_nodes.svg'
set encoding iso_8859_1
set xlabel 'Sockets'
set ylabel 'Efficiency (%)'
set xrange [4:32]
set xtics 4,4,32   # set an increment of 2
set yrange [0:100]
set ytics 0,10,100   # set an increment of 2
set nokey
set style data line
first(x) = ($0 > 0 ? base : base = x)
plot "weak_scaling_nodes.out" u 1:(first($2), base/$2*100.0)
