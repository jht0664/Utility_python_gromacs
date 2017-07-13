#!/bin/bash
# rename traj.*.xtc by reducing number by 1

INIT=$1
str1=".xtc"
str11="traj."

rm $str11$INIT$str1
let INIT=INIT+1


while [ -f "$str11$INIT$str1" ]
do 
	let TARGET=INIT-1
	str1trj=$str11$INIT$str1
	str2trj=$str11$TARGET$str1	
	mv $str1trj $str2trj
	let INIT=INIT+1
done 

