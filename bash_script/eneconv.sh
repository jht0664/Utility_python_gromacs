#!/bin/bash
# cat traj files in npt files located in ../#npt/
# $1 : initial number for folder
# $2 : final number for folder
str1="npt"
path="../"
space=" "
init=$1
final=$2
traj=""
filename="/ener.edr"

until [ $init -gt $final ]
do
	temp=$space$path$init$str1$filename
	traj=$traj$temp
	let init=init+1
done

gmx_mpi eneconv -f $traj -settime -o ener.edr
