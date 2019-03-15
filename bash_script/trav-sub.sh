#!/bin/bash
# traverse npt-subdirectories in a current directory
# for folder "#str"
# $1 : initial number
# $2 : final number
# $3 : string folder
# $4 : subfolder name
# $5, 6 ,7 : command lines in each folder


INIT=$1
FINAL=$2
str1=$3
subfolder=$4
command3="gmx_mpi energy"

until [ $INIT -gt $FINAL ]
do 
	strdir=$INIT$str1	
	cd $strdir
	cd $subfolder
	pwd
# do a command-line
	$5 $6 $7
#	echo $command2 | $command3
	let INIT=INIT+1
	cd ..
	cd ..
done 

