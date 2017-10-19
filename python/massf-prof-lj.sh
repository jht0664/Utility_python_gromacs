#!/bin/bash
# make trajectory using index file 
# $1 : nbins

a_select="../../../ndx/a.select"
b_select="../../../ndx/b.select"
mass_info="../../../ndx/select.mass"
mass_info2="../../../ndx/select.mass2"

echo $a_select
if [ -f "$a_select" ]
then 
	echo "$a_select found."
else
	echo "$a_select not found."
	exit 1
fi
echo $b_select
if [ -f "$b_select" ]
then
        echo "$b_select found."
else
        echo "$b_select not found."
        exit 1
fi
echo $mass_info
if [ -f "$mass_info" ]
then
        echo "$mass_info found."
else
        echo "$mass_info not found."
        exit 1
fi

## extract trajectory
## for bin size
#source activate mdaenv
python ~/Utility/python/massf-prof.py -i traj.trr -s topol.tpr -m $mass_info -select1 $a_select -select2 $b_select -nbin $1 -axis 2 -o a | tee a.massf.log 
python ~/Utility/python/massf-prof.py -i traj.trr -s topol.tpr -m $mass_info2 -select1 $b_select -select2 $a_select -nbin $1 -axis 2 -o b | tee b.massf.log
#source deactivate
