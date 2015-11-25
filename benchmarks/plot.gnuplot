set terminal svg enhanced size 500,500
set output 'file.svg'
set encoding iso_8859_1
set title 'This is the title of the graph'
set xlabel 'Generation'
set ylabel 'Average'
set nokey
set yrange [0:100]
plot 'stdevPlayer1.dat' using 1:($2*100) with lines, 'stdevPlayer2.dat' using 1:($3*100) with lines