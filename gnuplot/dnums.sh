#!/bin/bash
# cat traj files in npt files located in ../#npt/
# $1 : box-z length, nm
# $2 : bin size, nm
# $3 : total time, ns
# you should set followings:
init1="b"
init2="peo"
text1="[BF_4]^-"
text2="PEO(n=9)"
text3="pre-separation"
natom1="1"
natom2="30"
temp="430"
minmax1="2.5" # abs(cbr) for b
minmax2="0.8" # abs(cbr) for peo
blocklength="20" # block length (frame)
frametime="0.04" # ns

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

for ifile in $file11 $file12 $file21 $file22
do
	cp ~/Utility/gnuplot/dnums.plt $ifile
done
for ifile in $file31 $file32
do
	cp ~/Utility/gnuplot/dnums.acfs.plt $ifile
done

sed -i "s/ARG1/$input11/g" $file11
sed -i "s/ARG1/$input12/g" $file12
sed -i "s/ARG1/$input21/g" $file21
sed -i "s/ARG1/$input22/g" $file22
sed -i "s/ARG1/$input31/g" $file31
sed -i "s/ARG1/$input32/g" $file32

for ifile in $file11 $file12 $file21 $file22 $file31 $file32
do
	sed -i "s/BOXZZZ/$1/g" $ifile
	sed -i "s/BINZZZ/$2/g" $ifile
	sed -i "s/TOTALTIME/$3/g" $ifile
	sed -i "s/OUTPUT/$ifile$exe2/g" $ifile
	sed -i "s/TEMPERATURE/$temp/g" $ifile
	sed -i "s/INITCOORD/$text3/g" $ifile
	sed -i "s/BLOCKLL/$blocklength/g" $ifile
	sed -i "s/TIMETT/$frametime/g" $ifile
done

# for b
for ifile in $file11 $file21 $file31
do
	sed -i "s/NATOM/$natom1/g" $ifile
	sed -i "s/RANGE/$minmax1/g" $ifile
	sed -i "s/MOLECULE/$text1/g" $ifile
done
# for peo
for ifile in $file12 $file22 $file32
do
	sed -i "s/NATOM/$natom2/g" $ifile
	sed -i "s/RANGE/$minmax2/g" $ifile
	sed -i "s/MOLECULE/$text2/g" $ifile
done

# run gnuplot
command1="gnuplot "
$command1 $file11
$command1 $file12
$command1 $file21
$command1 $file22
$command1 $file31
$command1 $file32

