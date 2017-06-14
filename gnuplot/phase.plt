reset

# postscript terminal
set terminal postscript enhanced color \
    font 'Verdana,22' linewidth 2
set output 'OUTPUT'

#set xlabel "{Box z-position, nm}"
#set ylabel "Time, ns"

#unset key

# border
set style line 11 lc rgb '#808080' lt 2
set border 15 front ls 11
set tics nomirror out scale 1.25

set xl "Fraction"
set yl "Density"
set xtics 0,0.2,1 
set mxtics 5
set grid

plot 'phasediagram.gnuplot' i 0 u 1:2:3 w xe lw 3 lt 2 lc 1 title "Ref", 'phasediagram.gnuplot' i 1 u 1:2:3 w xe lw 3 lt 2 lc 1 notitle ,'phasediagram.gnuplot' i 4 u 1:2 w p lw 3 lt 5 lc 4 title "MD (T=1)", 'phasediagram.gnuplot' i 5 u 1:2 w p lw 3 lt 5 lc 4 notitle, 'phasediagram.gnuplot' i 2 u 1:2 w p lw 3 lt 5 lc 3 title "MD (T=2)", 'phasediagram.gnuplot' i 3 u 1:2 w p lw 3 lt 5 lc 3 notitle


#set xrange [0:boxz]
#set yrange [0:endtime]


