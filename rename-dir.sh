#!/bin/bash
# rename directory

INIT=$1
FINAL=$2
str1=$3
str2=$4

until [ $INIT -gt $FINAL ]
do 
	str1dir=$INIT$str1
	str2dir=$INIT$str2	
	mv $str1dir $str2dir
	let INIT=INIT+1
done >trav-npt.log

