#!/bin/bash
# run analysis for lj phase seprataion
# $1 : folder name
# $2 : new folder name

mkdir $2
cd $2
cp ../$1/confout.gro conf.gro
cp ../$1/grompp.mdp ./
cp ../$1/table* ./
cp ../$1/index.ndx ./
cp ../$1/run_eq.sh ./
cp ../$1/topol.top ./

sed -i "s/continuation             = no/continuation             = yes/g" grompp.mdp
sed -i "s/gen-vel                  = yes/gen-vel                  = no/g" grompp.mdp
sed -i "s/$1/$2/g" run_eq.sh

qsub run_eq.sh

