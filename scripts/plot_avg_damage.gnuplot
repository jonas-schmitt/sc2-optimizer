set terminal svg enhanced size 500,500
set output 'avg_damage.svg'
set encoding iso_8859_1
set xlabel 'Generation'
set ylabel 'Damage'
set yrange [0:100]
set ytics 0,10,100
set key bottom left
set style data line
plot "./results/avgPlayer1.dat" u 1:($2*100.0) title "Player 1" lc 3, "./results/avgPlayer2.dat" u 1:($2*100.0) title "Player 2" lc 1
