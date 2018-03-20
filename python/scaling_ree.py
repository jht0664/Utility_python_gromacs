#!/usr/bin/env python3
# ver 0.1 - copy from rdf_itf.py (v0.1) and modify codes on 2/3/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculation scaling of Ree of single chain')
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
	help='cut-off checking distance between atoms in a molecule (d_cutoff < d_neighbor_atoms: stop)')
parser.add_argument('-b', '--begin', default=-1, nargs='?', type=int,
	help='begining frame (-1: last half trajectory)')
parser.add_argument('-o', '--output', default='pol.ree.scal', nargs='?', 
	help='output filename for scaling of Ree file')
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
from scipy.spatial.distance import pdist

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
atomtxt = open(args.select).read()
#hjung.polymer.check_traj_connectivity(u,str(atomtxt),args.nmol,args.cutoff,'simple')

## data setting
select_mol = u.select_atoms(str(atomtxt))
if len(select_mol)%args.nmol != 0:
	raise ValueError("wrong # molecules, (args.nmol, select_mol) {} {} ".format(args.nmol, len(select_mol)))
n_deg = int(len(select_mol)/args.nmol)
print("assume all molecules has {} atoms".format(n_deg))
data_ree = np.zeros((args.nmol,n_deg-1))
#data_ree_vec = np.zeros((args.nmol,n_deg,3))

# make a list, the indices which are in the same lag
def list_idx_decr(start,end,init_step):
	#print("in {} {} {}".format(start,end,init_step))
	list_i = []
	while start < end:
		if start < 0:
			print(" something wrong (bug?)!")
			break
		list_i.append(int(start))
		init_step = init_step - 1
		start = start + init_step
	return list_i

# make list set (list of list of indices) where the elements (in fact, list of indices) are grouped by lag
# pair_data_size   number of pairs between any two points excluding duplicates 
#                  which is the same as the length of result array in scipy.spatial.distance.pdist
# n_data_points    number of datas you used for scipy.spatial.distance.pdist
#                  which is the same as the length of argument in scipy.spatial.distance.pdist
def main_list_idx(pair_data_size,n_data_points):
	# check validity of arguments
	#print(" main_list_idx:")
	expect_size = int((n_data_points-1)*n_data_points/2)
	if int(pair_data_size) != expect_size:
		raise ValueError(" Your arugments are wrong because {}(input) != {}(expect) based on {} ".format(pair_data_size,expect_size,n_data_points))
	list_set = []
	i_end = pair_data_size
	max_lag = n_data_points
	for i_start in range(n_data_points):
		#print(" lag: {}".format(i_start))
		i_end = i_end - i_start
		if i_end < i_start:
			break
		list_set.append(list_idx_decr(i_start,i_end,max_lag))
	return list_set

## read trajectory
i_frame = 0
imod = hjung.time.process_init()
list_sets = main_list_idx(int((n_deg-1)*n_deg/2),n_deg)
#print(list_sets)
print(len(list_sets))
for ts in u.trajectory[skip_frames:]:
	for i_mol in range(args.nmol):
		pair_dist = pdist(select_mol.positions[i_mol*n_deg:(i_mol+1)*n_deg],metric='euclidean')
		#print(pair_dist[list_sets[len(list_sets)-1]])
		#print(np.mean(pair_dist[list_sets[len(list_sets)-1]]))
		for i_lag in range(len(list_sets)):
			data_ree[i_mol,i_lag] = data_ree[i_mol,i_lag] + np.mean(pair_dist[list_sets[i_lag]])
	i_frame = i_frame + 1
	imod = hjung.time.process_print(i_frame, n_frames, imod)
#print(float(n_frames))
norm_data_ree = data_ree/float(n_frames)
norm_data_ree = np.transpose(norm_data_ree)

# save raw rg data file
np.savetxt(args.output, norm_data_ree, 
	header='time-averaged scaling of Ree of {} single chains'.format(args.nmol), fmt='%e', comments='# ')
np.save(args.output, norm_data_ree)

# save avg file
if args.nmol > 1:
	data_ree_avg = np.column_stack((np.mean(norm_data_ree, axis=1),np.std(norm_data_ree, axis=1)))
	np.savetxt(args.output+'.avg', data_ree_avg, 
		header='molecule-averaged scaling of Ree of {} chains'.format(args.nmol), fmt='%e', comments='# ')
	print(" saved average Ree files")

## timer
hjung.time.end_print(start_proc, start_prof)
