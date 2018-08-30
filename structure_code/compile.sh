#!/bin/bash

OMPINCLUDE=/usr/mpi/intel/openmpi-1.2.8/lib64/openmpi/
MKLLIB=/opt/intel/mkl/10.0.3.020/lib/em64t/

#OPT="-openmp -static -check bounds -check uninit -check format -warn declarations -traceback"  #-warn unused 

OPT="-openmp -static"

ifort $OPT -c glob_v.f90 routines.f90 pme.f90 structure_factor.f90 main_structure.f90  -I$MKLLIB -I$OMPINCLUDE

ifort $OPT *.o -L$MKLLIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$OMPINCLUDE -libompitv -o main_structure