#!/bin/bash
# make trajectory using index file 
# $1 : nbins
# $2 : tolerance

a_select="../../ndx/a.select"
b_select="../../ndx/b.select"

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

## extract trajectory
## for bin size
#~/anaconda3/bin/python3 ~/Utility/python/dens.py -select $peo_select -bin $1 -onum peo.nums -obin peo.bins -o peo.prob
#~/anaconda3/bin/python3 ~/Utility/python/dens.py -select $b_select -bin $1 -onum b.nums -obin b.bins -o b.prob
python3 ~/Utility/python/dn-auto.py -i traj.trr -s topol.tpr -select $b_select -nbin $1 -tol $2 -axis 2 -o b.dnums
python3 ~/Utility/python/dn-auto.py -i traj.trr -s topol.tpr -select $a_select -nbin $1 -tol $2 -axis 2 -o a.dnums


