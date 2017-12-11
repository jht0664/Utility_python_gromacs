#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 02/25/2017
# ver 0.2 - support pdb and dcd files for openmm on 5/8/2017
# ver 0.3 - support xtc trajectory files for Monte Carlo using "reduce_unitcells_3d_to_1d" on 6/6/2017
# ver 0.4 - support block average module on 6/26/2017 and remove out again.
# ver 0.5 - add some function: printing arguments and reduce lines as for default setting
# ver 1.0 - as for a template, in alignment, use step function, instead of using directly autocorrelation function.
#			This way can reflect the case that a mole fraction profile has a large flat region.
#			Previous way is sensitive for a little of sharp tip in a large flat region. on 11/28/2017
# ver 1.1 - separate making slab geometry and alignment functions on 11/28/2017
# ver 1.2 - save additional file for iframe when two more layers make (this is an exceptional case)
#			and remove multilayer trajectories on 11/30/2017
# ver 1.3 - save files for domain sizes and mole fraction and total density in center of domain on 12/6/2017
# ver 1.4 - copy from conv_align.py on 12/10/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='1D alignment using shift array')
## args
parser.add_argument('-i', '--input', default='traj.mass', nargs='?', 
	help='raw 1D profile (npy file format) and exclude .npy in argument for alignment')
parser.add_argument('-s', '--shift', default='traj.mass.conv', nargs='?', 
	help='shift array (npy file format) and exclude .npy in argument')
parser.add_argument('-del', '--delete', default='traj.mass.removeframe', nargs='?',
	help='index of frames to exclude')
parser.add_argument('-o', '--output', default='.align', nargs='?', 
	help='output surfix for aligned profiles')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.4')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np
import matplotlib
matplotlib.use('Agg') # avoid to show figures when running bash shell script
import matplotlib.pyplot as plt
from scipy import stats
import copy

# default for args
args.output = args.input + args.output # save aligned mass fraction profiles
args.input = args.input + '.npy'
args.shift = args.shift + '.npy'
args.delete = args.delete + '.npy'

## timer
start_proc, start_prof = hjung.time.init()

## load data files
mass_1d_t = np.load(args.input)
shift = np.load(args.shift)
removeframe = np.load(args.delete)
if len(mass_1d_t) != len(shift):
	raise ValueError("the size of data files are different.")
nbin = len(mass_1d_t[0])

## alignment using shift array and remove using removeframe arrat
for iframe in range(len(mass_1d_t)):
	shift_array = mass_1d_t[iframe]
	mass_1d_t[iframe] = np.roll(shift_array, shift[iframe])
align_mass_1d_t = np.delete(mass_1d_t,removeframe,axis=0)	

## write
np.savetxt(args.output, align_mass_1d_t, 
	header='%d, %d, aligned mass by shift arrary' \
	%(len(align_mass_1d_t),nbin), fmt='%f', comments='# ')
#np.save(args.omassf, align_massfrac_1d_t) 

## timer
hjung.time.end_print(start_proc, start_prof)