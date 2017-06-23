#!/bin/bash
# copy important files from ../#npt/
# $1 : initial number for folder
# $2 : final number for folder
# $3 : nbin
# $4 : tolerance or block length
str1="npt"
path="../"
space=" "
init=$1
final=$2
traj=""

cp $path$init$str1"/conf.gro" ./
cp $path$final$str1"/confout.gro" ./
cp $path$final$str1"/grompp.mdp" ./
cp $path$init$str1"/topol.tpr" ./

~/Utility/eneconv.sh $init $final
~/Utility/energy.sh
~/Utility/trjcat.sh $init $final

# for Gromacs
~/Utility/python/dn-auto.sh $3 $4 | tee dn-auto.log
~/Utility/python/massf-prof.sh $3 | tee massf.log
# for OpenMM
#~/Utility/python/dn-auto-openmm.sh $3 $4 | tee dn-auto.log
#~/Utility/python/massf-prof-openmm.sh $3 | tee massf.log

# dnum avg

# massf avg 
#python ~/Utility/python/savetxt-avg.py -i peo.massf.align -b 7500 -tol $4 | tee massf-avg.log

#gnuplot ./massfit.log
