#!/usr/bin/env python3
# ver 0.1 - copy from massf_1d.py (v1.1) and modify codes on 1/22/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='partial rdf(s) at interfaces using various domain info. from domain.py')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr structure file')
parser.add_argument('-select1', '--select1', nargs='?', 
	help='a file1 with a command-line for select_atoms in MDAnalysis')
parser.add_argument('-select2', '--select2', nargs='?', 
	help='a file2 with a command-line for select_atoms in MDAnalysis')
parser.add_argument('-nbin', '--nbin', nargs='?', type=int,
	help='number of bins when you did conv. alignmnet')
parser.add_argument('-d', '--domain', default='a.massf.domain', nargs='?',
	help='.npz file from domain.py')
parser.add_argument('-di', '--domain_i', default='a.massf.domain.dic', nargs='?',
	help='input text file for determining domain.py')
parser.add_argument('-hnbin', '--hist_nbin', default=20, nargs='?', type=int,
	help='number of bins for rdfs')
parser.add_argument('-hmax', '--hist_max', default=-1.0, nargs='?', type=float,
	help='the maximum distance for rdf (if negative, use half box_x)')
parser.add_argument('-temp', '--temp', default=20, nargs='?', type=int,
	help='generate tmp file to save intermediate data')
parser.add_argument('-o', '--output', default='.rdf', nargs='?', 
	help='output prefix for rdf files')
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
import numpy as np
import math
import multiprocessing as mp
#from multiprocessing import Pool
#import os
#import resource
import time
import copy # copy array
import ast  # reading a file of dictionary list
from scipy.spatial.distance import cdist

# default for args
odomain = args.domain + args.output
idomain = args.domain + '.npy'

## timer
start_proc, start_prof = hjung.time.init()

#if args.n_proc < 0:
#	args.n_proc = mp.cpu_count() - 1
#if args.n_proc > mp.cpu_count():
#	raise ValueError(" args.n_proc is greater than number of cpus in this node")

## load domain data files
domain_info = np.load(idomain)
n_frames = len(domain_info)
#if n_frames != len(domain_info):
#	raise ValueError(" Wrong trajectory length. Check if your trajectory file is changed.")
iframes_list = domain_info[:,0]
align_shift = domain_info[:,2]
step_up_down = np.zeros((n_frames,2))
step_up_down[:,0] = domain_info[:,3]
step_up_down[:,1] = domain_info[:,4]
domain_size = domain_info[:,5]

def reduce_traj_output(file_gro, file_xtc, file_select1, file_select2, list_frames):
	coordinates1, coordinates2, unit_cells = hjung.io.read_trr_3d_select2(file_gro, file_xtc, file_select1, file_select2, 'pos')
	unit_cells = hjung.array.convert_unitcell_3d(unit_cells, file_gro, file_xtc)
	return coordinates1[list_frames], coordinates2[list_frames], unit_cells[list_frames]

# read a topology and a trajectory using module MDAnalysis with selection
coordinates1, coordinates2, unit_cells = reduce_traj_output(args.structure, args.input, args.select1, args.select2, iframes_list)

# read domain_dic
with open(args.domain_i,'r') as file_dic:
	dic_domain = ast.literal_eval(file_dic.read())
	print(" Here is dic_domain: {}".format(dic_domain))
n_rdfs = int(np.amax(np.array(list(dic_domain.values())))) + 1

# interface position calculation
itf = np.zeros(2) # [0]: (down -> up) massf interface, [1]: (up -> down) massf interface
itf[0] = int(args.nbin/4.0)
itf[1] = int(args.nbin*3.0/4.0)

def select_atoms_indices_1d(coord,r_min,r_max):
	right_indice_atoms = np.where(coord > r_min)
	left_indice_atoms = np.where(coord < r_max)
	return np.intersect1d(right_indice_atoms,left_indice_atoms)

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

# save rdf files
def save_rdf(rdf,n_rdfs,bin_edges,count_frames,odomain,prefix):
	rdf_tmp = copy.copy(rdf)
	# normalization against volume of bins
	for i_bin in range(len(bin_edges)-1):
		vol_bin = (4.0/3.0)*math.pi*((bin_edges[i_bin+1])**3-(bin_edges[i_bin])**3)
		## normalization against number of frames
		for i_rdf in range(n_rdfs):
			rdf_tmp[i_rdf][i_bin] = rdf[i_rdf][i_bin]/vol_bin/float(count_frames[i_rdf])

	# bin position (center)
	hist_nbins = int(len(bin_edges)-1)
	hist_bin_pos = np.zeros(hist_nbins,dtype=float)
	for ibin in range(hist_nbins):
		hist_bin_pos[ibin] = (bin_edges[ibin+1] + bin_edges[ibin])/2.0
	# save number histogram trajectory
	for i_rdf in range(n_rdfs):
		dataset = np.column_stack((hist_bin_pos,rdf_tmp[i_rdf]))
		np.savetxt(odomain+str(prefix)+str(i_rdf), dataset, 
			header='radial distribution function at domain median +- {} bins with {} frames'.format(i_rdf,count_frames), fmt='%e', comments='# ')
		np.save(odomain+str(prefix)+str(i_rdf), dataset)
	print(" saved (temp) rdf files")

# calculation rdfs
count_frames = np.zeros(n_rdfs, dtype=np.int)
hist_nbins = args.hist_nbin
rdf = np.zeros((n_rdfs,hist_nbins), dtype=np.float)
current_i = 0
box_z = unit_cells[0][2] # assume unit cell dimension is fixed during simulation (NVT)
bin_size = box_z/float(args.nbin)
box_x_half = unit_cells[0][0]/2.0
print("bin size on z = {}".format(bin_size))
print("bin_x = {}".format(unit_cells[0][0]))
if args.hist_max < 0.0:
	args.hist_max = box_x_half
elif args.hist_max > box_x_half:
	raise ValueError(" wrong args.hist_max")
print("set hist_max = {}".format(args.hist_max))

## multiprocessing but so slow
#def calc_dist(atom1,atom2):
#	global coordinates1_tmp
#	global coordinates2_tmp
#	return np.linalg.norm(coordinates1_tmp[atom1] - coordinates2_tmp[atom2])
#
#def pool_calc_dist_hist(nproc, pair_list, nbins, max_dist):
#	pool = Pool(nproc)
#	t = time.time()
#	distances = pool.starmap(calc_dist, pair_list)
#	#print(" total {} pairs".format(len(dist_list)))
#	print(" time = {}".format(time.time()-t))
#	hist_dist, bin_edges = np.histogram(distances,bins=nbins,range=(0.0,max_dist),density=False)
#	return hist_dist, bin_edges

for iframe in range(n_frames):
	i_rdf = int(dic_domain[domain_size[iframe]])
	pos_itf = ((itf - align_shift[iframe])%args.nbin)*bin_size
	for up_down in range(2):
		new_coord_1_z = coordinates1[iframe][:,2] - pos_itf[up_down] + box_z/2.0 # to set interface at center of box
		new_coord_2_z = coordinates2[iframe][:,2] - pos_itf[up_down] + box_z/2.0 # to set interface at center of box
		new_coord_1_z = hjung.coord.pbc_nojump_1d(new_coord_1_z, unit_cells[iframe][2])
		new_coord_2_z = hjung.coord.pbc_nojump_1d(new_coord_2_z, unit_cells[iframe][2])
		# selection selection1 of atoms within a range
		#print("{} {}".format(step_up_down[iframe][up_down], itf[up_down]))
		t_itf = abs(step_up_down[iframe][up_down] - itf[up_down])*bin_size/2.0 # (half_thickness between interface and acf 50% point) / 2
		#print(" selection1 range = [{},{}] = {} +- {}".format(box_z/2.0 - t_itf, box_z/2.0 + t_itf,box_z/2.0, t_itf))
		select_atoms_1 = select_atoms_indices_1d(new_coord_1_z, box_z/2.0 - t_itf, box_z/2.0 + t_itf)
		print(" {} frame: selected {} atoms in select1 ".format(iframe,len(select_atoms_1)))
		# selection atoms pool of selection2 atoms
		select_atoms_2 = select_atoms_indices_1d(new_coord_2_z, box_z/2.0 - t_itf - box_x_half, box_z/2.0 + t_itf + box_x_half)
		#print(" selection2 range = [{},{}]".format(box_z/2.0 - t_itf - box_x_half,box_z/2.0 + t_itf + box_x_half))
		print(" {} frame: selected {} atoms in select2".format(iframe,len(select_atoms_2)))
		
		## rdf of A-B near interface using multiprocessing (slower by 400 times than scipy)
		#comb_list = [(i,j) for i in select_atoms_1 for j in select_atoms_2]
		#hist_dist, bin_edges = pool_calc_dist_hist(args.n_proc, comb_list,hist_nbins,args.hist_max)
		# use cdist
		#t = time.time()
		result = cdist(coordinates1[iframe][select_atoms_1], coordinates2[iframe][select_atoms_2], 'euclidean')
		dist_tmp = np.ndarray.flatten(result)
		dist_tmp_reduced = dist_tmp[np.nonzero(dist_tmp)]
		#print(" time = {}".format(time.time()-t))
		hist_dist, bin_edges = np.histogram(dist_tmp_reduced,bins=hist_nbins,range=(0.0,args.hist_max),density=False)
		# save rdfs
		rdf[i_rdf] = rdf[i_rdf] + hist_dist/len(select_atoms_1) # normalizaed for atom1
		count_frames[i_rdf] = count_frames[i_rdf] + 1
	if (current_i != 0) and (current_i%10 == 0):
		save_rdf(rdf,n_rdfs,bin_edges,count_frames,odomain,'.tmp.')
	
save_rdf(rdf,n_rdfs,bin_edges,count_frames,odomain,'.')

## timer
hjung.time.end_print(start_proc, start_prof)