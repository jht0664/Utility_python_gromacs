#!/bin/bash
# delete every second line in a file
# $1 : initial number for folder
# $2 : final number for folder
# $3 : n-th line you want to take
rm b.nums
rm peo.nums

str1="npt"
init=$1
final=$2

until [ $init -gt $final ]
do
	folder=$init$str1
	awk 'NR%'$3'==0' ../$folder/b.nums >> b.nums
	awk 'NR%'$3'==0' ../$folder/peo.nums >> peo.nums
	let init=init+1
done


