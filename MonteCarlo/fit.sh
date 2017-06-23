#/bin/bash
# $1 = nbins

ls traj.*.xtc > wc.log
lframe=$(wc wc.log | awk '{print $1}')
let lfile=lframe-1
../../trjcat-mc.sh 0 $lfile
rm wc.log
cp ../../select* ./
python3 ~/Utility/python/massf-prof.py -i traj_out.xtc -s conf.gro -m select.mass -select1 select.txt -select2 select2.txt -nbin $1 -axis 2 -o a.massf > massf.log
btime=$(grep 'frames' massf.log | awk '{ print $4}')
let btime=btime/2
python3 ~/Utility/python/savetxt-avg.py -i a.massf.align -tol 0 -b $btime
cp ../../fit.plot ./
sed -i "s/NBINS/$1/g" fit.plot
gnuplot fit.plot 
cat fit.log
