set terminal svg enhanced size 500,500
set output 'processors_runtime.svg'
set encoding iso_8859_1
set xlabel 'Processors'
set ylabel 'Runtime (s)'
set xrange [1:8]
set yrange [0:1200]
set xtics 1,1,8
set nokey
set style data line
plot "weak_scaling_processors.out" u 1:($2/1000.0)
