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
	help='selection of end1 of each molecule')
parser.add_argument('-select2', '--select2', nargs='?',
	help='selection of end2 of each molecule')
parser.add_argument('-nmol', '--nmol', nargs='?', type=int,
	help='# molecules')
parser.add_argument('-cutoff', '--cutoff', default=0.0, nargs='?', type=float,
	help='cut-off checking distance between atoms in a molecule (d_cutoff < d_neighbor_atoms: stop)')
parser.add_argument('-b', '--begin', default=-1, nargs='?', type=int,
	help='begining frame (-1: last half trajectory)')
parser.add_argument('-o', '--output', default='pol', nargs='?', 
	help='output prefix filename for Ree files (.ree)')
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

args.output = args.output + '.ree'

## timer
start_proc, start_prof = hjung.time.init()

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
n_frames = n_frames - skip_frames
#print("goal n_frames = {}".format(n_frames))
atomtxt1 = open(args.select1).read()
atomtxt2 = open(args.select2).read()
print(" your selection1: {}".format(atomtxt1))
print(" your selection2: {}".format(atomtxt2))
#hjung.polymer.check_traj_connectivity(u,str(atomtxt),args.nmol,args.cutoff,'simple')

## data setting
data_ree = np.zeros((n_frames,args.nmol))
data_ree_vec = np.zeros((n_frames,args.nmol,3))
select_mol1 = u.select_atoms(str(atomtxt1))
select_mol2 = u.select_atoms(str(atomtxt2))
if len(select_mol1)%args.nmol != 0:
	raise ValueError("wrong # molecules, (args.nmol, select_mol) {} {} ".format(args.nmol, len(select_mol1)))
if len(select_mol2)%args.nmol != 0:
	raise ValueError("wrong # molecules, (args.nmol, select_mol) {} {} ".format(args.nmol, len(select_mol2)))

## read trajectory
i_frame = 0
imod = hjung.time.process_init()
for ts in u.trajectory[skip_frames:]:
	data_ree[i_frame] = distances.dist(select_mol1,select_mol2)[2]
	#vec_tmp = select_mol2.positions - select_mol1.positions
	#print(np.sqrt(np.sum(np.dot(vec_tmp[0],vec_tmp[0]))))
	#print(data_ree[i_frame][0])
	data_ree_vec[i_frame] = select_mol2.positions - select_mol1.positions
	i_frame = i_frame + 1
	imod = hjung.time.process_print(i_frame, n_frames, imod)
print("read total {} frames".format(i_frame))

# save raw rg data file
np.savetxt(args.output, data_ree, 
	header='Ree (mean = {} +- {} with {} frames'.format(np.mean(data_ree),np.std(data_ree),n_frames), fmt='%f', comments='# ')
np.save(args.output, data_ree)
print("Ree = {:.3f} +- {:.3f}".format(np.mean(data_ree),np.std(data_ree)))
print(" saved ree files")

# ree.vec
np.save(args.output+str('.vec'), data_ree_vec.flatten())

# save avg file
data_ree_tavg = np.column_stack((np.mean(data_ree, axis=0),np.std(data_ree, axis=0)))
np.savetxt(args.output+'.tavg', data_ree_tavg, 
	header='averaged Ree for each molecule with {} frames'.format(n_frames), fmt='%f', comments='# ')
print(" saved average Ree files")

## timer
hjung.time.end_print(start_proc, start_prof)