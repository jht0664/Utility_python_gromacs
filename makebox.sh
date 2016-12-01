#!/bin/bash
# insert bmim 

COUNTER=0
echo COUNTER $COUNTER
until [ $COUNTER -gt 20 ]; do
	let SEEDNUM=COUNTER+4412
	let COUNTERT=COUNTER+1
	gmx_mpi insert-molecules -f re$COUNTER.gro -ci ../../../FF/peo_9.gro -nmol 40 -try 1000 -seed $SEEDNUM -o re$COUNTERT.gro -allpair
	let COUNTER=COUNTER+1
done

gmx_mpi insert-molecules -f re5.gro -ci ../../../../FF/peo_9.gro -nmol 40 -try 2000 -allpair -o re6.gro -seed 5334


