#!/bin/bash
# generate folders like "1str", "2str", ...
# $1 : initial number for folder
# $2 : final number for folder
# $3 : folder name
str1=$3
init=$1
final=$2

until [ $init -gt $final ]
do
	let old=init-1
	oldfolder=$old$str1
	newfolder=$init$str1
	let old2=init-2
	old2folder=$old2$str1
	let new2=init+1
	new2folder=$new2$str1
	mkdir $newfolder
	#cp $oldfolder/conf.gro $newfolder/
	#cp $oldfolder/run_eq.sh $newfolder/
	#cp $oldfolder/table* $newfolder/
	#cp 1temp/grompp* $newfolder/
	#cp $oldfolder/index.ndx $newfolder/
	#cp $oldfolder/topol.top $newfolder/
        #cd $newfolder/
	#value=$( echo "120.2717*(1.6+0.2*$init)" | bc )
	#sed -i "s/216.48906/$value/g" grompp.mdp
	#sed -i "s/$newfolder/$new2folder/g" run_eq.sh
	#sed -i "s/$oldfolder/$newfolder/g" run_eq.sh
	#sed -i "s/$old2folder/$oldfolder/g" run_eq.sh
	#qsub run_eq.sh
	#cd ../
	let init=init+1
done


