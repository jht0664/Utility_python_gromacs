#!/bin/bash -f

qstat -u htjung | grep htjung > qstat.imsi
test=$(cat qstat.imsi | cut -c1-7)
echo
echo "**********************************"
echo "  qdel all jobs $test"
echo "**********************************"
echo

qdel $test
rm qstat.imsi

