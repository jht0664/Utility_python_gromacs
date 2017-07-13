#!/bin/bash
# cat traj files in npt files located in ../#npt/
# $1 : initial number for folder
# $2 : final number for folder
path="./"
space=" "
init=$1
final=$2
traj=""
filename1="/traj.trr"
filename2="traj."
filename3=".xtc"

until [ $init -gt $final ]
do
	temp=$space$path$filename2$init$filename3
	traj=$traj$temp
	let init=init+1
done

gmx trjcat -f $traj -settime -o traj_out.xtc

