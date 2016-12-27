#!/bin/bash
# traverse npt-subdirectories in a current directory
# with command-line $1 (ex. ~/Utility/n-prof.sh) 

INIT=8
FINAL=19
str1="npt"

until [ $INIT -gt $FINAL ]
do 
	strdir=$INIT$str1	
	cd $strdir
	pwd
# do a command-line
	$1 $2
	let INIT=INIT+1
	cd ..
done >trav-npt.log

