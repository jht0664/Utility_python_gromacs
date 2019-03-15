#!/bin/bash
# clean all intermediate files
# if we use -all option, delete conf.gro 

if [ $1 == "-all" ]
then
	rm conf.gro
fi

rm \#*
rm confout.gro 
rm ener.edr energy.xvg
rm msd.xvg
rm log md.log mdout.mdp topol.tpr
rm *.out
rm *.cpt
rm traj.*
rm *.pdb

echo "All clean!"

