#!/bin/bash
# make trajectory using index file 
# $1 : nbins

peo_select="../../../ndx/peo.select"
b_select="../../../ndx/b.select"
mass_info="../../../ndx/select.mass"
mass_info2="../../../ndx/select2.mass"

echo $peo_select
if [ -f "$peo_select" ]
then 
	echo "$peo_select found."
else
	echo "$peo_select not found."
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
python3 ~/Utility/python/massf-prof.py -i traj.trr -s topol.tpr -m $mass_info -select1 $peo_select -select2 $b_select -nbin $1 -axis 2 -o peo.massf
python3 ~/Utility/python/massf-prof.py -i traj.trr -s topol.tpr -m $mass_info2 -select1 $b_select -select2 $peo_select -nbin $1 -axis 2 -o b.massf 

