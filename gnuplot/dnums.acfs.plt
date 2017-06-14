#!/usr/bin/gnuplot
##
## Plotting a color map using the default Matlab palette
##
## AUTHOR: Hagen Wierstorf
# ARG1 : filename to plot
reset

#boxx = 4.2 #nm
#boxy = 4.2 #nm
boxz = BOXZZZ #nm
binsize = BINZZZ #nm
endtime = TOTALTIME #ns
blocklength = BLOCKLL #frames
frametime = TIMETT #ns
#npart = 1 # 1 for bf4 and 30 for pep-9mer

abscbr = 1
set cbrange[-abscbr:abscbr]
set title 'Density fluctuation of MOLECULE at TEMPERATURE K starting INITCOORD'

# postscript terminal
set terminal postscript enhanced color \
    font 'Verdana,18' linewidth 2
set output 'OUTPUT'

set xlabel "{Distance, nm}"
set ylabel "{Time, ns}"
set cblabel "{Normalized Density Fluctuation}"

unset key

# border
set style line 11 lc rgb '#808080' lt 1
set border 15 front ls 11
set tics nomirror out scale 0.75

set xrange [-boxz/2:boxz/2]
set yrange [0:endtime]

#disable colorbar tics
set cbtics scale 0
set cbtics 1

set pm3d map
load '/home/hjung52/Utility/gnuplot/palette.plt' 
splot ARG1 u ($1*binsize-boxz/2):($2*blocklength*frametime):3 matrix

