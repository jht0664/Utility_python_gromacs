#!/bin/bash
# make trajectory using index file 

ndxfile="../../../ndx/"$1".ndx"
echo $1
if [ $1 != "ll5l" ]; then
	if [ "$1" != "ll6l" ]; then
		if [ "$1" != "ll8l" ] ;then
			if [ "$1" != "small" ] ;then
				echo "Use proper argument 1"
			fi
		fi
	fi
fi

echo $ndxfile

if [ -f "$ndxfile" ]
then 
	echo "$ndxfile found."
else
	echo "$ndxfile not found."
	exit 1
fi

# extract trajectory (5 6 for OPLS,small(OPLS-400k) or 7 8 for ll5l-peo10,peo40 or 8 7 for ll6l-peo10)
echo 5 | gmx_mpi traj -n $ndxfile -oxt peo.xtc
echo 6 | gmx_mpi traj -n $ndxfile -oxt b.xtc


# do density profile program
program="./D-prof.x"
if [ -f "$program" ]
then
        echo "$program found."
else
        echo "$program not found."
        cp ~/Utility/density/compile-linux/D-Prof.x ./
fi

./D-Prof.x -pf peo -slice 0.4 -block 1
./D-Prof.x -pf b -slice 0.4 -block 1

# plot peo.Z.dens.avg and b.Z.dens.avg
