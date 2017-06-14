reset 

set title 'Mass fraction of PEO at different force fields'

# postscript terminal
set terminal postscript enhanced color \
    font 'Verdana,18' linewidth 2
set output 'mass.avg.eps'

set xlabel "{z/box-z}"
set ylabel "Mass Fraction of PEO"

# border
set style line 11 lc rgb '#808080' lt 1
set border 15 front ls 11
set tics nomirror out scale 0.75

set xrange [0:1]
set yrange [0:1.3]

f(x) = wp+0.5*(wr-wp)*(tanh((x-b+c)/lamda)-tanh((x-b-c)/lamda))
wr=0.996585;wp=0.121754;b=0.503065;c=0.235427;lamda=0.037778

g(x)=wp1+0.5*(wr1-wp1)*(tanh((x-b1+c1)/lamda1)-tanh((x-b1-c1)/lamda1))
wr1=0.952305;wp1=0.149675;b1=0.502845;c1=0.230386;lamda1=0.0556816

plot 'Research2/phase-bf4/ll5l/480K/separate/data/peo.massf.align.avg' u ($0/70):1:2 w ye lt 4 lc 6 title "OPLS data at 480 K", f(x) lw 4 lc 6 title "OPLS data fit", 'ua_il_peo_lcst_test/600K/peo.massf.align.avg' u ($0/90):1:2 w ye lt 7 lc 7 lw 1 title "CG data at 600 K", g(x) lw 4 lc 7 title "CG data fit"


