#!/bin/bash
# traverse npt-subdirectories in a current directory
# with command-line $1 (ex. ~/Utility/n-prof.sh) 

INIT=$1
FINAL=$2
str1="npt"

until [ $INIT -gt $FINAL ]
do 
	strdir=$INIT$str1	
	cd $strdir
	pwd
# do a command-line
	$3 $4 $5
	let INIT=INIT+1
	cd ..
done >trav-npt.log

