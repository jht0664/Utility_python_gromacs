#!/bin/bash
# check validity of trajectory file 
#  to avoid magic number error and zero trajectories

ls traj.*.xtc > wc.log
lframe=$(wc wc.log | awk '{print $1}')
let lfile=lframe-1
rm wc.log

start1=0
until [ $start1 -gt $lfile ]
do
	echo "... Checking ${start1}th file ..."
	../../trjcat-mc.sh $start1 $start1 >& check-trj.imsi &
	sleep 1
	error1=$(grep "Magic Number Error" check-trj.imsi | wc | awk '{print $1}')
	error2=$(grep "Inconsistent results" check-trj.imsi | wc | awk '{print $1}')
	let errort=error1+error2

	if [ $errort -ne 0 ]
	then
		echo "remove ${start1}-th file"
		../../rename-trj.sh $start1
		let lfile=lfile-1
	else
		let start1=start1+1
	fi
done
rm check-trj.imsi

echo "Done!"

