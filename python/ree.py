#!/usr/bin/env python3
# ver 0.1 - copy from rdf_itf.py (v0.1) and modify codes on 2/3/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculation end-to-end distance of molecules you select')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr structure file')
parser.add_argument('-select1', '--select1', nargs='?',
	help='selection1 of end atom of each molecule')
parser.add_argument('-select2', '--select2', nargs='?',
	help='selection atom2')
parser.add_argument('-nmol', '--nmol', nargs='?', type=int,
	help='# molecules')
parser.add_argument('-b', '--begin', default=-1, nargs='?', type=int,
	help='begining frame (-1: last half trajectory)')
parser.add_argument('-o', '--output', default='pol.ree', nargs='?', 
	help='output filename for Ree files')
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
from MDAnalysis.analysis import distances
import numpy as np
from scipy.spatial.distance import euclidean

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

atom1txt = open(args.select1).read()
atom2txt = open(args.select2).read()

u = mda.Universe(args.structure,args.input)
n_frames = len(u.trajectory)
if args.begin == -1:
	skip_frames = int(n_frames/2)
	print(" skip {} frames".format(skip_frames))
else:
	skip_frames = args.begin
if args.begin >= n_frames:
	raise ValueError("wrong args.begin because of > n_frames")
data_ree = np.zeros((n_frames,args.nmol))
i_frame = 0
imod = process_init()
end1 = u.select_atoms(str(atom1txt))
end2 = u.select_atoms(str(atom2txt))
if (args.nmol != len(end1)) or (args.nmol != len(end2)):
	raise ValueError("wrong # molecules, (args.nmol, end1, end2) {} {} {} ".format(args.nmol, len(end1), len(end2)))
origin = np.zeros(3)
for ts in u.trajectory:
	diff_pos = end1.positions - end2.positions
	pos_wrap = hjung.coord.pbc_nojump_3d_test(diff_pos, ts._unitcell)
	i_atom = 0
	for i_atom in range(args.nmol):
		#print(euclidean(i_pos, origin))
		data_ree[i_frame,i_atom] = euclidean(pos_wrap[i_atom], origin)
		i_atom = i_atom + 1
	i_frame = i_frame + 1
	imod = process_print(i_frame, n_frames, imod)
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
np.savetxt(args.output, data_ree, 
	header='Rg (mean = {} +- {} with {} frames'.format(np.mean(data_ree[skip_frames:]),np.std(data_ree[skip_frames:]),n_frames-skip_frames), fmt='%f', comments='# ')
np.save(args.output, data_ree)
print("Ree = {} +- {}".format(np.mean(data_ree[skip_frames:]),np.std(data_ree[skip_frames:])))
print(" saved ree files")

# save avg file
data_ree_tavg = np.column_stack((np.mean(data_ree, axis=0),np.std(data_ree, axis=0)))
np.savetxt(args.output+'.tavg', data_ree_tavg, 
	header='averaged Ree for each molecule with {} frames'.format(n_frames), fmt='%f', comments='# ')
data_ree_mavg = np.column_stack((np.mean(data_ree, axis=1),np.std(data_ree, axis=1)))
np.savetxt(args.output+'.mavg', data_ree_mavg, 
	header='averaged Ree for each frame with {} molecules'.format(args.nmol), fmt='%f', comments='# ')
print(" saved average Ree files")

## timer
hjung.time.end_print(start_proc, start_prof)