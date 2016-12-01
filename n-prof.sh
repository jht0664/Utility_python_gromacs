#!/bin/bash
# make trajectory using index file 

echo "bin size = $1"
if [ $# -ne 1 ]; then
	echo "$0: usage: n-prof.sh bin_size(float)"
	echo "we assume bin_size is 2.0 A"
	$1=2.0
fi

peo_select="../../../ndx/peo.select"
b_select="../../../ndx/b.select"

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

# extract trajectory
~/anaconda3/bin/python3 ~/Utility/python/dens.py -select $peo_select -bin $1 -onum peo.nums -obin peo.bins -o peo.prob

~/anaconda3/bin/python3 ~/Utility/python/dens.py -select $b_select -bin $1 -onum b.nums -obin b.bins -o b.prob

