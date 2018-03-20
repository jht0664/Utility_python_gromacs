#!/usr/bin/env python3
# ver 0.1 - copy from rdf_itf.py (v0.1) and modify codes on 2/3/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculation isotropic form vector')
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
parser.add_argument('-o', '--output', default='pol.form_iso', nargs='?', 
	help='output filename for form factor (iso) files')
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
from scipy.spatial.distance import pdist
import math


## timer
start_proc, start_prof = hjung.time.init()


## read files
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

## set range of q
box_min = np.amin(u.trajectory[0]._unitcell) # assume we use NVT ensemble simulation
min_q = 1.0/(box_min/2.0)
max_q = 1.0 # the unit for distance is A, then unit of q is A^(-1)
list_q = np.arange(min_q,max_q,min_q)
n_qs = len(list_q)
#print("list of q's")
#print(list_q)


## size of dist list
n_comb = math.factorial(n_deg)/math.factorial(n_deg-2)/2 # number of combination of distances, "nC2"


## distance calculations
data_form_iso = np.zeros((args.nmol,n_qs))
#data_form_iso_int = np.zeros((args.nmol,n_qs))
i_frame = 0
imod = hjung.time.process_init()
skip_steps = int(10)
for ts in u.trajectory[skip_frames::skip_steps]:
	for i_mol in range(args.nmol):
		# check the validity of atomic positions 
		start_iatom = i_mol*n_deg
		end_iatom = (i_mol+1)*n_deg
		dist = pdist(select_mol.positions[start_iatom:end_iatom])
		i_qs = 0
		for iq in list_q:
			data_form_iso[i_mol,i_qs] = data_form_iso[i_mol,i_qs] + np.sum(np.divide(np.sin(dist*iq),dist*iq))/float(n_comb)
			#data_form_iso_int[i_mol,i_qs] = data_form_iso_int[i_mol,i_qs] + (np.sum(np.divide(np.sin(dist*iq),dist*iq))/float(n_comb))**2
			i_qs += 1
	i_frame = i_frame + 1
	imod = hjung.time.process_print(i_frame, n_frames/skip_steps, imod)

data_form_iso = np.divide(data_form_iso,float(i_frame))
data_form_iso_avg = np.mean(data_form_iso,axis=0)
data_form_iso_std = np.std(data_form_iso,axis=0)


# transpose
data_form_iso = np.transpose(data_form_iso)
data_form_iso_avg = np.transpose(data_form_iso_avg)
data_form_iso_std = np.transpose(data_form_iso_std)
data_set = np.column_stack((list_q,data_form_iso))
data_set_avg = np.column_stack((list_q,data_form_iso_avg,data_form_iso_std))


## save
np.savetxt(args.output, data_set,
	header='qs, F(q)s - form factor - of {} molecules'.format(args.nmol), fmt='%e', comments='# ')
np.save(args.output, data_set)
np.savetxt(args.output+str('.avg'), data_set_avg,
	header='qs, avg F(q), std F(q) of all {} molecules'.format(args.nmol), fmt='%e', comments='# ')
np.save(args.output+str('.avg'), data_set_avg)


## timer
hjung.time.end_print(start_proc, start_prof)