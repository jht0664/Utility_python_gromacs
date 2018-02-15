#!/usr/bin/env python3
# ver 0.1 - copy from rdf_itf.py (v0.1) and modify codes on 2/3/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculation end-to-end distance of molecules you select using vector')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr structure file')
parser.add_argument('-select', '--select', nargs='?',
	help='selection of each molecule')
parser.add_argument('-nmol', '--nmol', nargs='?', type=int,
	help='# molecules')
parser.add_argument('-cutoff', '--cutoff', default=0.0, nargs='?', type=float,
	help='cut-off distance between atmos in a molecule (d_cutoff < d_neighbor_atoms: stop)')
parser.add_argument('-b', '--begin', default=-1, nargs='?', type=int,
	help='begining frame (-1: last half trajectory)')
parser.add_argument('-o', '--output', default='pol.ree_vec', nargs='?', 
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

## read files
u = mda.Universe(args.structure,args.input)
n_frames = len(u.trajectory)
if args.begin == -1:
	skip_frames = int(n_frames/2)
	print(" skip {} frames".format(skip_frames))
else:
	skip_frames = args.begin
if args.begin >= n_frames:
	raise ValueError("wrong args.begin because of > n_frames")
atomtxt = open(args.select).read()

## data setting
data_ree = np.zeros((n_frames,args.nmol))
dist_max = 0.0
dist_min = 10.0
select_mol = u.select_atoms(str(atomtxt))
if len(select_mol)%args.nmol != 0:
	raise ValueError("wrong # molecules, (args.nmol, select_mol) {} {} ".format(args.nmol, len(select_mol)))
n_deg = int(len(select_mol)/args.nmol)
print("assume all molecules has {} atoms".format(n_deg))

## read trajectory
i_frame = 0
imod = process_init()
origin = np.zeros(3)
for ts in u.trajectory:
	for i_mol in range(args.nmol):
		# check the validity of atomic positions to calculate Ree
		for i_atom in range(n_deg-1):
			dist = euclidean(select_mol.positions[i_mol*n_deg+i_atom], select_mol.positions[i_mol*n_deg+i_atom+1])
			if dist > args.cutoff:
				print("maybe due to the wrapped trajectory or too small cutoff.")
				raise RuntimeError("pos[{}}][{}] = {} and {}, but dist = {} > cutoff {}".format(
					i_mol,i_atom,select_mol.positions[i_atom], select_mol.positions[i_atom], dist, args.cutoff))
			if dist > dist_max:
				dist_max = dist
			if dist < dist_min:
				dist_min = dist
		data_ree[i_frame,i_mol] = euclidean(select_mol.positions[i_mol*n_deg], select_mol.positions[(i_mol+1)*n_deg-1])
	i_frame = i_frame + 1
	imod = process_print(i_frame, n_frames, imod)

print("Ree ranges [{},{}]".format(dist_min, dist_max))

# save raw rg data file
np.savetxt(args.output, data_ree, 
	header='Rg (mean = {} +- {} with {} frames'.format(np.mean(data_ree),np.std(data_ree),n_frames), fmt='%f', comments='# ')
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