#!/bin/bash
# make .plt file for mass fraction
# $1 : box-x and box-y (they should be same)
# $2 : box-z length, nm
# $3 : nbin
# $4 : nframe
# $5 : total time, ns
# $6 : temperature
# you should set followings:
init1="peo" # input file name prefix
text1="PEO(n=9)" # y lable for molecules
text3="PEO(n=9) mass fraction" # title for plot

echo $1 : box-x and box-y nm
echo $2 : avg box-z length, nm
echo $3 : nbin
echo $4 : nframe
echo $5 : total time, ns
echo $6 : temperature

# preset
file1=".massf"
file2=".massf.align"
exe1=".plt"
exe2=".eps"
comm="'"
input11=$comm$init1$file1$comm #b.massf
input21=$comm$init1$file2$comm #b.massf.align

for ifile in $file11 $file12 
do
	cp ~/Utility/gnuplot/massf.plt $ifile
done

sed -i "s/ARG1/$input11/g" $file11
sed -i "s/ARG1/$input12/g" $file12

for ifile in $file11 $file12 
do
	sed -i "s/BOXXXX/$1/g" $ifile	
	sed -i "s/BOXYYY/$1/g" $ifile
	sed -i "s/BOXZZZ/$2/g" $ifile
	sed -i "s/NBINNN/$3/g" $ifile
	sed -i "s/OUTPUT/$ifile$exe2/g" $ifile
	sed -i "s/CUSTOMTITLE/$text3/g" $ifile
	sed -i "s/TIMETT/$5/g" $ifile
	sed -i "s/FRAMETT/$4/g" $ifile
	sed -i "s/TEMPERATURE/$6/g" $ifile
	sed -i "s/MOLECULE/$text1/g" $ifile
done

# run gnuplot
command1="gnuplot "
$command1 $file11
$command1 $file12

