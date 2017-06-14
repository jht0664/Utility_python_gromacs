#!/usr/bin/gnuplot
##
## Plotting a color map using the default Matlab palette
##
## AUTHOR: Hagen Wierstorf
# ARG1 : filename to plot
reset

#boxx = 8 #nm
#boxy = 8 #nm
#boxz = BOXZZZ #nm 
#binsize = BINZZZ #nm
#endtime = TOTALTIME #ns
#blocklength = BLOCKLL #frames
#frametime = TIMETT #ns
#npart = NATOM # 1 for bf4 and 30 for pep-9mer

#abscbr = RANGE # 2.5 for BF4, 1 for PEO
set title 'Density Profile of PEO-9 at Different Temperatures'

# postscript terminal
set terminal postscript enhanced color \
    font 'Verdana,18' linewidth 2
set output 'peo.dens.eps'

set xlabel "{x/box-z}"
set ylabel "{Number density, /nm^3}"

set key

# border
set style line 11 lc rgb '#808080' lt 1
set border 3 front ls 11
set tics nomirror out scale 0.75

set xrange [0:1]
set yrange [0:1.4]

#disable colorbar tics
plot 'll5l/600K/separate/2npt/peo.dnums.align.avg' u ($0/120):($1/30/0.2044875/4.2/4.2):($2/30/0.2044875/4.2/4.2) w ye lt 1 title "600 K", \
'll5l/600K/separate/2npt/peo.dnums.align.avg' u ($0/120):($1/30/0.2044875/4.2/4.2):($2/30/0.2044875/4.2/4.2) w l lw 2 lt 1 notitle, \
'll5l/540K/separate/data/peo.dnums.align.avg' u ($0/114):($1/30/0.2019894/4.2/4.2):($2/30/0.2019894/4.2/4.2) w ye lt 2 title "540 K", \
'll5l/540K/separate/data/peo.dnums.align.avg' u ($0/114):($1/30/0.2019894/4.2/4.2):($2/30/0.2019894/4.2/4.2) w l lw 2 lt 2 notitle, \
'll5l/480K/separate/data/peo.dnums.align.avg' u ($0/110):($1/30/0.1975785/4.2/4.2):($2/30/0.1975785/4.2/4.2) w ye lt 3 title "480 K", \
'll5l/480K/separate/data/peo.dnums.align.avg' u ($0/110):($1/30/0.1975785/4.2/4.2):($2/30/0.1975785/4.2/4.2) w l lw 2 lt 3 notitle, \
'll5l/440K/separate/data/peo.dnums.align.avg' u ($0/70):($1/30/0.2983445/4.2/4.2):($2/30/0.2983445/4.2/4.2) w ye lt 4 title "440 K", \
'll5l/440K/separate/data/peo.dnums.align.avg' u ($0/70):($1/30/0.2983445/4.2/4.2):($2/30/0.2983445/4.2/4.2) w l lw 2 lt 4 notitle, \
'll2l/420k/separate/data/peo.dnums.align.avg' u ($0/56):($1/30/0.1470802/4/4):($2/30/0.1470802/4/4) w ye lt 5 title "420 K", \
'll2l/420k/separate/data/peo.dnums.align.avg' u ($0/56):($1/30/0.1470802/4/4):($2/30/0.1470802/4/4) w l lw 2 lt 5 notitle
