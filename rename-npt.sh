#!/bin/bash
# rename directory

INIT=$1
FINAL=$2
str1="NPAT"
str2="npt"

until [ $INIT -gt $FINAL ]
do 
	str1dir=$str1$INIT
	str2dir=$INIT$str2	
	mv $str1dir $str2dir
	let INIT=INIT+1
done >trav-npt.log

