#!/bin/bash

BBCOUNTER=1111
upperlim=40

for ((i=26; i<=upperlim; i++)); do
   echo The counter is "$i"
   let CCOUNTER=i+1
   /home/htjung/gromacs5/bin/gmx_mpi_d insert-molecules -f re$i.gro -ci ../scale-ref/no3.gro -nmol 200 -try 10 -seed $BBCOUNTER -o re$CCCOUNTER.gro
   let BBCOUNTER=BBCOUNTER+1
done



