#!/bin/bash
# cat traj files in npt files located in ../#npt/
# $1 : initial number for folder
# $2 : final number for folder
# $3 : initial filename
str1="npt"
init=$1
final=$2
traj=""
left="["
right="]"
hyphen="-"
filename=$left$init$hyphen$final$right
space=" "

until [ $init -gt $final ]
do
	temp=$space$init$str1
	folder=$folder$temp
	let init=init+1
done

tar -cvf $3$filename.tar $folder

