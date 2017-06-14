#!/usr/bin/gnuplot
##
## Plotting a color map using the default Matlab palette
##
## AUTHOR: Hagen Wierstorf
# ARG1 : filename to plot
reset

boxx = 4.2 #nm
boxy = 4.2 #nm
boxz = BOXZZZ #nm 
binsize = BINZZZ #nm
endtime = TOTALTIME #ns
blocklength = BLOCKLL #frames
frametime = TIMETT #ns
npart = NATOM # 1 for bf4 and 30 for pep-9mer

abscbr = RANGE # 2.5 for BF4, 1 for PEO
set cbrange[-abscbr:abscbr]
set title 'Number Density of MOLECULE at TEMPERATURE K starting INITCOORD'

# postscript terminal
set terminal postscript enhanced color \
    font 'Verdana,18' linewidth 2
set output 'OUTPUT'

set xlabel "{Box z-position, nm}"
set ylabel "Time, ns"
set cblabel "{Number density, /nm^3}"

unset key

# border
set style line 11 lc rgb '#808080' lt 1
set border 15 front ls 11
set tics nomirror out scale 0.75

set xrange [0:boxz]
set yrange [0:endtime]

#disable colorbar tics
set cbtics scale 0
set cbtics RANGE

set pm3d map
load '/home/hjung52/Utility/gnuplot/palette.plt' 
splot ARG1 u ($1*binsize):($2*blocklength*frametime):($3/npart/binsize/boxx/boxy) matrix

