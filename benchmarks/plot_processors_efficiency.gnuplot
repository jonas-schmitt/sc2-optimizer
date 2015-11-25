set terminal svg enhanced size 500,500
set output 'processors_efficiency.svg'
set encoding iso_8859_1
set xlabel 'Processors'
set ylabel 'Efficiency (%)'
set xrange [1:8]
set xtics 1,1,8
set yrange [0:100]
set ytics 0,10,100
set nokey
set style data line
first(x) = ($0 > 0 ? base : base = x)
plot "weak_scaling_processors.out" u 1:(first($2), base/$2*100.0)
