#!/usr/bin/env python3
# ver 0.1 - copy from rdf_itf.py (v0.1) and modify codes on 3/28/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculation center of mass of polymers you select')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr structure file')
parser.add_argument('-select', '--select', nargs='?',
	help='selection of each molecule')
parser.add_argument('-nmol', '--nmol', nargs='?', type=int,
	help='# molecules')
parser.add_argument('-cutoff', '--cutoff', default=15.0, nargs='?', type=float,
	help='cut-off checking distance between atoms in a molecule (d_cutoff < d_neighbor_atoms: stop)')
parser.add_argument('-b', '--begin', default=-1, nargs='?', type=int,
	help='begining frame (-1: last half trajectory)')
parser.add_argument('-unit', '--unitcell', default='NO', nargs='?',
	help='save unit cell info? (YES/NO)')
parser.add_argument('-o', '--output', default='pol', nargs='?', 
	help='output prefix filename for COM files (.com or like pol.com) ')
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

args.output = args.output + '.com'

## timer
start_proc, start_prof = hjung.time.init()

## read files
u = mda.Universe(args.structure,args.input)
n_frames = len(u.trajectory)
#print(n_frames)
if args.begin == -1:
	skip_frames = int(n_frames/2)
	print(" skip {} frames".format(skip_frames))
else:
	skip_frames = args.begin
if args.begin >= n_frames:
	raise ValueError("wrong args.begin because of > n_frames")
n_frames = n_frames - skip_frames
atomtxt = open(args.select).read()
#if ('.trr' in args.input) or ('.xtc' in args.input):
#	print(" you seem to use Gromacs trajectory, thus we are going to check connectivity within intramolecular atoms.")
#	hjung.polymer.check_traj_connectivity(u,str(atomtxt),args.nmol,args.cutoff,'simple')
print("#"*15)
print(" Be aware that this program does not consider periodic boundary condition.")
print(" Please check wrapping option like gmx trjconv -pbc mol for Gromacs trajectory.")
print("#"*15)

## data setting
data_com = np.zeros((n_frames,args.nmol,3))
data_unitcell = np.zeros((n_frames,3))
select_mol = u.select_atoms(str(atomtxt))
if len(select_mol)%args.nmol != 0:
	raise ValueError("wrong # molecules, (args.nmol, select_mol) {} {} ".format(args.nmol, len(select_mol)))
n_deg = int(len(select_mol)/args.nmol)
print("assume all molecules has {} atoms".format(n_deg))

## read trajectory
from numpy import linalg as LA
i_frame = 0
imod = hjung.time.process_init()
for ts in u.trajectory[skip_frames:]:
	for i_mol in range(args.nmol):
		com_tmp = select_mol[i_mol*n_deg:(i_mol+1)*n_deg].center_of_mass() # But sometimes com is beyond box dimension
		data_com[i_frame,i_mol,0] = com_tmp[0] - u.dimensions[0]*int(round(com_tmp[0]/u.dimensions[0]-0.50))
		data_com[i_frame,i_mol,1] = com_tmp[1] - u.dimensions[1]*int(round(com_tmp[1]/u.dimensions[1]-0.50))
		data_com[i_frame,i_mol,2] = com_tmp[2] - u.dimensions[2]*int(round(com_tmp[2]/u.dimensions[2]-0.50))
		#print(i_frame,i_mol)
		#if (i_mol == 921):
		#	print(com_tmp)
		#	print(data_com[i_frame,i_mol])
		#	print(u.dimensions)
		#	#print(select_mol[i_mol*n_deg:(i_mol+1)*n_deg].position)
		## check connectivity
		#for ii in range(i_mol*n_deg,(i_mol+1)*n_deg-1):
		#	if LA.norm(select_mol[ii+1].position-select_mol[ii].position) > 20:
		#		print("wrong? {} {} {} {}".format(LA.norm(select_mol[ii+1].position-select_mol[ii].position),i_mol,ii-i_mol*n_deg,i_frame))
		if (data_com[i_frame,i_mol,2] > u.dimensions[2]) or (data_com[i_frame,i_mol,2] < 0):
			print("somthing strange {} {} {} {}".format(i_frame, i_mol, i_mol*n_deg, data_com[i_frame,i_mol,2]))
	data_unitcell[i_frame] = u.dimensions[0:3]
	i_frame = i_frame + 1
	imod = hjung.time.process_print(i_frame, n_frames, imod)
#if data_unitcell[i_frame-1][0] < 1.0:
#	print("last frame has null trajectory, so reduce nframes")
#	i_frame = i_frame - 1
#data_com = data_com[0:i_frame]
#data_unitcell = data_unitcell[0:i_frame]
print("read total {} frames".format(i_frame))

#### TEST: calc. com histograms
#unit_cells_1d = data_unitcell[:,2]
#com_1d = data_com[:,:,2]
#com_hist_1d_t, bin_1d_t = hjung.analyze.histo_t_1d_nbin(com_1d, unit_cells_1d, 142) 
#np.savetxt('test.com.histo', com_hist_1d_t, fmt='%f')

# save com data file
flat_data = data_com.flatten()
np.savetxt(args.output, flat_data, 
	header='COM (mean = {} +- {} with {} frames and {} molecules'.format(np.mean(data_com),np.std(data_com),n_frames,args.nmol), fmt='%f', comments='# ')
np.save(args.output, flat_data)
if 'YES' in args.unitcell:
	file_unitcell = 'unit_cell'
	np.savetxt(file_unitcell, data_unitcell, 
		header='XYZ box dimensions with {} frames '.format(i_frame), fmt='%f', comments='# ')
	np.save(file_unitcell, data_unitcell)
	print("saved unitcell data files")

## timer
hjung.time.end_print(start_proc, start_prof)