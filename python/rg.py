#!/usr/bin/env python3
# ver 0.1 - copy from rdf_itf.py (v0.1) and modify codes on 2/3/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculation radius of gyration (Rg) of molecules you select')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr structure file')
parser.add_argument('-b_rid', '--b_rid', nargs='?', type=int,
	help='beginning of residu id of target molecules (starting 1)')
parser.add_argument('-nres', '--nres', nargs='?', type=int,
	help='number of residues of a single molecule')
parser.add_argument('-nmol', '--nmol', nargs='?', type=int,
	help='number of molecules you select')
parser.add_argument('-o', '--output', default='pol.rg', nargs='?', 
	help='output filename for Rg files')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import sys
sys.path.append('/home/htjung/Utility/python/')
import hjung
from hjung import *
import MDAnalysis as mda
import numpy as np

## timer
start_proc, start_prof = hjung.time.init()

# initialize processing bar
# output: process tiem and wall time
# Example: mod_frame = process_init()
def process_init():
	print("io.process_init: ")
	return 10 # initial mode number = 10

# print process bar
# input: itime, current time
#        ftime, final time
#        mod_time, frequency of printing.
# output: mod_time, new or old printing frequency 
# Example: mod_frame = process_print(itime+1, ftime, mod_frame)
def process_print(itime, ftime, mod_time):
	if itime%mod_time == 0:
		print("... {0} th frame reading ({1:.0%}) ...".format(itime,itime/ftime))
		if (itime/mod_time)%10 == 0:
			mod_time = mod_time*10
	return mod_time

u = mda.Universe(args.structure,args.input)
n_frames = len(u.trajectory)
data_rg = np.zeros((n_frames,args.nmol))
i_frame = 0
imod = process_init()
for ts in u.trajectory:
	i_rid = args.b_rid
	for i_mol in range(args.nmol):
		mol = u.select_atoms('resid '+str(i_rid)+'-'+str(i_rid+args.nres-1))
		data_rg[i_frame,i_mol] =  mol.radius_of_gyration()
		i_rid = i_rid + args.nres
	i_frame = i_frame + 1
	imod = process_print(i_frame, n_frames, imod)
print("read from resid {} to resid {} for {} frames".format(args.b_rid,i_rid-1,n_frames))
#
#      xt = x(i) - x(j)
#      yt = y(i) - y(j)
#      zt = z(i) - z(j)
#      XT = XT - box(1)*DNINT(XT/box(1))
#      YT = YT - box(2)*DNINT(YT/box(2))
#      ZT = ZT - box(3)*DNINT(ZT/box(3))
#      RT = (XT*XT+YT*YT+ZT*ZT) !*pos_scaling2
#      DO k=1, n_dv  ! count number of overlapping pairs (finally, double counting)
#        IF(RT*scale_length(k) < DCOMP) THEN
#          ompv(id,k) = ompv(id,k) + 1 
#          !write(*,*) id,k,ompv(id,k)
#        ENDIF
#      ENDDO

# save raw rg data file
np.savetxt(args.output, data_rg, 
	header='Rg (mean = {} +- {} with {} frames'.format(np.mean(data_rg),np.std(data_rg),n_frames), fmt='%f', comments='# ')
np.save(args.output, data_rg)
print("Rg = {} +- {}".format(np.mean(data_rg),np.std(data_rg)))
print(" saved rg files")

# save avg file
data_rg_tavg = np.column_stack((np.mean(data_rg, axis=0),np.std(data_rg, axis=0)))
np.savetxt(args.output+'.tavg', data_rg_tavg, 
	header='averaged Rg for each molecule with {} frames'.format(n_frames), fmt='%f', comments='# ')
data_rg_mavg = np.column_stack((np.mean(data_rg, axis=1),np.std(data_rg, axis=1)))
np.savetxt(args.output+'.mavg', data_rg_mavg, 
	header='averaged Rg for each frame with {} molecules'.format(args.nmol), fmt='%f', comments='# ')
print(" saved average rg files")

## timer
hjung.time.end_print(start_proc, start_prof)