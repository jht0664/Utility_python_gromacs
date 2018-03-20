#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculation radius of gyration using MDAnalysis')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr structure file')
parser.add_argument('-select', '--select', nargs='?',
	help='selection of each molecule')
parser.add_argument('-nmol', '--nmol', nargs='?', type=int,
	help='# molecules')
parser.add_argument('-b', '--begin', default=-1, nargs='?', type=int,
	help='begining frame (-1: last half trajectory)')
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

# read trajectory
u = mda.Universe(args.structure,args.input)
n_frames = len(u.trajectory)
skip_frames = 0
if args.begin == -1:
	skip_frames = int(n_frames/2)
	print(" skip {} frames".format(skip_frames))
else:
	skip_frames = args.begin
if args.begin >= n_frames:
	raise ValueError("wrong args.begin because of > n_frames")
n_frames = n_frames - skip_frames
atomtxt = open(args.select).read()
#hjung.polymer.check_traj_connectivity(u,str(atomtxt),args.nmol,1.8,'random')
select_mol = u.select_atoms(str(atomtxt))
if len(select_mol)%args.nmol != 0:
	raise ValueError("wrong # molecules, (args.nmol, select_mol) {} {} ".format(args.nmol, len(select_mol)))
n_deg = int(len(select_mol)/args.nmol)
print("assume {} atoms you select per molecule".format(n_deg))

# calculation of Rg
data_rg = np.zeros((n_frames,args.nmol))
i_frame = 0
imod = hjung.time.process_init()
for ts in u.trajectory[skip_frames:]:
	for i_mol in range(args.nmol):
		mol = select_mol.atoms[n_deg*i_mol:n_deg*(i_mol+1)]
		data_rg[i_frame,i_mol] =  mol.radius_of_gyration()
	i_frame = i_frame + 1
	imod = hjung.time.process_print(i_frame, n_frames, imod)

# save raw rg data file
np.savetxt(args.output, data_rg, 
	header='Rg for each molecules (mean = {} +- {}) with {} frames'.format(np.mean(data_rg),np.std(data_rg),n_frames), fmt='%f', comments='# ')
np.save(args.output, data_rg)
print("average Rg = {} +- {}".format(np.mean(data_rg),np.std(data_rg)))

# save avg file
data_rg_tavg = np.column_stack((np.mean(data_rg, axis=0),np.std(data_rg, axis=0)))
np.savetxt(args.output+'.avg', data_rg_tavg, 
	header='averaged Rg for each molecule with {} frames'.format(n_frames), fmt='%f', comments='# ')


## timer
hjung.time.end_print(start_proc, start_prof)