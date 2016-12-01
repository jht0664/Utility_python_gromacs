#!/bin/bash
# make trajectory using index file 

ndxfile="../../../ndx/"$1".ndx"
echo $1
if [ $1 != "ll5l" ]; then
	if [ "$1" != "ll6l" ]; then
		if [ "$1" != "ll8l" ] ;then
			echo "Use proper argument 1"
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

# extract trajectory
echo 7 | gmx_mpi traj -n $ndxfile -oxt peo.xtc
echo 8 | gmx_mpi traj -n $ndxfile -oxt b.xtc

# do density profile program
program="D-prof.x"
if [ -f "$program" ]
then
        echo "$program found."
else
        echo "$program not found."
        cp ~/Utility/density/compile-linux/D-Prof.x ./
fi

./D-Prof.x -pf peo -slice 0.3 -block 10
./D-Prof.x -pf b -slice 0.3 -block 10

# plot peo.Z.dens.avg and b.Z.dens.avg
