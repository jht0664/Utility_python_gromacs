#!/bin/bash
# make traj 

ndxfile_for_peo="../peo.ndx"
ndxfile_for_b="../b.ndx"
if [ -f "$ndxfile_for_peo" ]
then 
	echo "$ndxfile_for_peo found."
else
	echo "$ndxfile_for_peo not found."
	exit 1
fi
if [ -f "$ndxfile_for_b" ]
then
        echo "$ndxfile_for_b found."
else
        echo "$ndxfile_for_b not found."
        exit 1
fi

# extract trajectory
echo 5 | gmx_mpi traj -n $ndxfile_for_peo -oxt peo.xtc
echo 5 | gmx_mpi traj -n $ndxfile_for_b -oxt b.xtc

# do density profile program
program="D-prof.x"
if [ -f "$program" ]
then
        echo "$program found."
else
        echo "$program not found."
        cp ~/Utility/code/gromacs.densityProfile/compile-linux/D-Prof.x ./
fi

./D-Prof.x -pf peo -slice 0.2 -block 10
./D-Prof.x -pf b -slice 0.2 -block 10

# plot peo.Z.dens.avg and b.Z.dens.avg
