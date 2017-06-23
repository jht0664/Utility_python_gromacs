#!/bin/bash
# run analysis for lj phase seprataion
# $1 : nbin
# $2 : tolerance or block length

## for Gromacs
~/Utility/python/dn-auto-lj.sh $1 | tee dn-auto.log
~/Utility/python/massf-prof-lj.sh $1 | tee massf.log

# massf avg 
btime=$(grep 'frames' massf.log | awk '{ print $4}')
#declare -i btime
let btimes=btime/2
python ~/Utility/python/savetxt-avg.py -i a.massf.align -b $btimes -tol $2 | tee a.massf-avg.log
cp ~/Utility/gnuplot/massfit.plot ./
sed -i "s/NBINS/$1/g" massfit.plot
gnuplot massfit.plot
cat fit.log

#python ~/Utility/python/savetxt-avg.py -i b.massf.align -b $1 -tol $2 | tee b.massf-avg.log

## dnum avg
#~/Utility/python/dnum-avg-lj.sh $1 $2
#cp ~/Utility/gnuplot/dnumsfit.plot ./dnumsfit-a.plot
#sed -i "s/NBINS/$1/g" dnumsfit-a.plot
#sed -i "s/COM1/a/g" dnumsfit-a.plot
#gnuplot dnumsfit-a.plot
#cat fit.log
#cp ~/Utility/gnuplot/dnumsfit.plot ./dnumsfit-b.plot
#sed -i "s/NBINS/$1/g" dnumsfit-b.plot
#sed -i "s/COM1/b/g" dnumsfit-b.plot
#gnuplot dnumsfit-b.plot
#cat fit.log

#gnuplot ./massfit.log
