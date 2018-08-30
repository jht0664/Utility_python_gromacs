# Axes style

set style data histogram
set style histogram cluster gap 1
set style fill solid border lt -1

set key right top vertical Right noreverse noenhanced autotitle nobox
set boxwidth 0.92 absolute

#set xtics border rotate by -45 autojustify
set xtics norangelimit

set style line 101 lc rgb '#000000' lt 1
set style line 102 lc rgb '#808080' lt 1

set border 31 front ls 101 lw 2 
#or set border 3 back ls 101

set ytics mirror in scale 1
set xtics mirror in scale 1
set ytics font ",26"
set xtics font ",26"
