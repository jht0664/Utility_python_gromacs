#!/bin/bash
# run analysis for lj phase seprataion
# $1 : nbin
# $2 : tolerance or block length

# for Gromacs
#~/Utility/python/dn-auto-lj.sh $1 | tee dn-auto.log
#~/Utility/python/massf-prof-lj.sh $1 | tee massf.log

# massf avg 
python ~/Utility/python/savetxt-avg.py -i a.massf.align -b $1 -tol $2 | tee a.massf-avg.log
python ~/Utility/python/savetxt-avg.py -i b.massf.align -b $1 -tol $2 | tee b.massf-avg.log

# dnum avg
~/Utility/python/dnum-avg-lj.sh $1 $2

#gnuplot ./massfit.log
