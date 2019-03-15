#!/bin/bash
# do all interactions using gmx rdf

cp ./../../../ndx/rdf.ndx ./
INIT=20000
FINAL=30000

for REF in HI CMI HWI HWI1 CUI CSI CTI PEOH
do
	for GROUP in F O
	do
		gmx_mpi rdf -f traj.trr -s topol.tpr -n rdf.ndx -o rdf.$REF.$GROUP.xvg -b $INIT -e $FINAL -ref $REF -sel $GROUP
	done
done


