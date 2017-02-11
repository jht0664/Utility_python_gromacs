#!/bin/bash
# generate rdf parallelization
# $1 : initial frame of trajectory
# $2 : final frame of trajectory
ncpus=16
title="int-rdf3d_"
init=$1
final=$2
filename="rdf3d.sh."

rm $filename* # remove pervious old child shell files
rm *.plot3d.*
rm slurm-*.out

let range=final-init+1
if [ $range%$ncpus == 0 ]
then
	let block=range/ncpus
else
	let block=range/ncpus+1
fi
echo "each process has the range, " $block

if [ $block -le $ncpus ]
then
	echo "the time block range is very short."
	echo "please decrease the ncpus"
	exit 0
fi


until [ $init -gt $final ]
do
	input=$filename$init
	cp ~/Utility/python/interaction_3d.sh $input
	title_num=$title$init
	sed -i "s/TITLE/$title_num/g" $input
	sed -i "/echo/i cd $(pwd)" $input # any better?	
	sed -i "s/INIT/$init/g" $input
	
	let fframe=init+block-1
	if [ $fframe -gt $final ]
	then
		let fframe=final
	fi
	sed -i "s/FINAL/$fframe/g"  $input
	qsub $input
	let init=init+block
done


