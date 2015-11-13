set terminal svg enhanced size 500,500
set output 'stdev_health.svg'
set encoding iso_8859_1
set xlabel 'Generation'
set ylabel 'Health'
set yrange [0:100]
set ytics 0,10,100
set key top left
set style data line
plot "./results/stdevPlayer1.dat" u 1:($3*100.0) title "Player 1" lc 3, "./results/stdevPlayer2.dat" u 1:($3*100.0) title "Player 2" lc 1
