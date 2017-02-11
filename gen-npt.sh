#!/bin/bash
# generate npt folder and queue
# $1 : initial number for folder
# $2 : final number for folder
str1="npt"
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
	cp $oldfolder/run_eq.sh $newfolder/
	cd $newfolder/
	sed -i "s/$newfolder/$new2folder/g" run_eq.sh
	sed -i "s/$oldfolder/$newfolder/g" run_eq.sh
	sed -i "s/$old2folder/$oldfolder/g" run_eq.sh
	cd ../
	let init=init+1
done


