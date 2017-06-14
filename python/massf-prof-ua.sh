#!/bin/bash
# make trajectory using index file 
# $1 : nbins

peo_select="../peo.select"
b_select="../b.select"
mass_info="../select.mass"

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
~/anaconda2/bin/python2 ~/Utility/python/massf-prof.py -i md_nvt_ns.dcd -s md_nvt_init.pdb -m $mass_info -select1 $peo_select -select2 $b_select -nbin $1 -axis 2 -o peo.massf


