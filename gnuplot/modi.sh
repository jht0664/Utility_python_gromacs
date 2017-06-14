#!/bin/bash
# cat traj files in npt files located in ../#npt/
# $1 : box-z length, nm
# $2 : bin size, nm
# $3 : total time, ns
# you should set followings:
init1="b"
init2="peo"

# preset
file1=".dnums"
file2=".dnums.align"
file3=".dnums.acfs"
exe1=".plt"
exe2=".eps"
comm="'"
input11=$comm$init1$file1$comm #b.dnums
input12=$comm$init2$file1$comm #peo.dnums
input21=$comm$init1$file2$comm #b.dnums.align
input22=$comm$init2$file2$comm #peo.dnums.align
input31=$comm$init1$file3$comm #b.dnums.acfs
input32=$comm$init2$file3$comm #peo.dnums.acfs
file11=$init1$file1$exe1 #b.dnums.plt
file12=$init2$file1$exe1 #peo.dnums.plt
file21=$init1$file2$exe1 #b.dnums.align.plt
file22=$init2$file2$exe1 #peo.dnums.align.plt
file31=$init1$file3$exe1 #b.dnums.acfs.plt
file32=$init2$file3$exe1 #peo.dnums.acfs.plt

for ifile in $file11 $file12 $file21 $file22 $file31 $file32
do
#	sed -i "s/boxx = 4.2/boxx = 4/g" $ifile
#	sed -i "s/boxy = 4.2/boxy = 4/g" $ifile
	sed -i "s/blocklength = 10/blocklength = 20/g" $ifile
done

#for ifile in $file11 $file21 
#do
#        sed -i "s/abscbr = /abscbr = 2.5 #/g" $ifile
#	sed -i "s/set cbtics /set cbtics 2.5 #/g" $ifile
#done

#for ifile in $file12 $file22
#do
#        sed -i "s/abscbr = /abscbr = 0.8 #/g" $ifile
#	sed -i "s/set cbtics /set cbtics 0.8 #/g" $ifile
#done

# run gnuplot
command1="gnuplot "
$command1 $file11
$command1 $file12
$command1 $file21
$command1 $file22
$command1 $file31
$command1 $file32

