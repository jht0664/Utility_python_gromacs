#!/bin/bash
# energies for NPAT ensemble 
index="1 2 3 4 5 6 7 8 11 19 21 22"
#index="1 2 3 4 5 6 7 9 12 21 23"
echo $index | gmx_mpi energy  >& energy.log

