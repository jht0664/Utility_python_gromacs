#!/usr/bin/gnuplot
##
## Plotting a color map using the default Matlab palette
##
## AUTHOR: Hagen Wierstorf

reset

set pm3d map
load '/home/hjung52/Utility/gnuplot/palette.plt'

set xlabel '{/Helvetica-Oblique box-z position, nm}'
set ylabel '{/Helvetica-Oblique time, ns}'

unset key

# border
set style line 11 lc rgb '#808080' lt 1
set border 3 front ls 11
set tics nomirror out scale 0.75

#disable colorbar tics
set cbtics scale 0


