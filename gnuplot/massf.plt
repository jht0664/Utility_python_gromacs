#!/usr/bin/gnuplot
##
## Plotting a color map using the default Matlab palette
##
## AUTHOR: Hagen Wierstorf
# ARG1 : filename to plot
reset
boxx = BOXXXX #nm
boxy = BOXYYY #nm
boxz = BOXZZZ #nm 
nbin = NBINNN #nm
frametime = FRAMETT #each frame time
endtime = TIMETT #ns

abscbr = 1 #
set cbrange[-abscbr:abscbr]
set title 'CUSTOMTITLE at TEMPERATURE'

# postscript terminal
set terminal postscript enhanced color \
    font 'Verdana,18' linewidth 2
set output 'OUTPUT'

set xlabel "{Box z-position, nm}"
set ylabel "Time, ns"
set cblabel "{MOLECULE mass fraction}"

unset key

# border
set style line 11 lc rgb '#808080' lt 1
set border 15 front ls 11
set tics nomirror out scale 0.75

set xrange [0:boxz]
set yrange [0:endtime]

#disable colorbar tics
set cbtics scale 0
set cbtics 1

set pm3d map
load '/home/hjung52/Utility/gnuplot/palette.plt' 
splot ARG1 u ($1/nbin):($2*frametime):($3) matrix

