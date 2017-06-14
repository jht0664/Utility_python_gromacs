#!/bin/bash
# traverse npt-subdirectories in a current directory
# for folder "#str"
# $1 : initial number
# $2 : final number
# $3 : string folder
# $4, 5 ,6 : command lines in each folder

INIT=$1
FINAL=$2
str1=$3
command3="gmx_mpi energy"

until [ $INIT -gt $FINAL ]
do 
	strdir=$INIT$str1	
	cd $strdir
	pwd
# do a command-line
	$4 $5 $6
#	echo $command2 | $command3
	let INIT=INIT+1
	cd ..
done >trav-npt.log

