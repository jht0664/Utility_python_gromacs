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
# ver 1.2 - copy from massf_1d.py on 12/10/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='1D mass profile using counting number of particles ')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr structure file')
parser.add_argument('-select', '--select', nargs='?', 
	help='a file with a command-line for select_atoms in MDAnalysis')
parser.add_argument('-nbin', '--nbin', nargs='?', type=int,
	help='number of bins')
parser.add_argument('-axis', '--axis', default=2, nargs='?', type=int,
	help='which axis for histogram (x axis (0), y axis (1), z axis (2))')
parser.add_argument('-o', '--output', default='traj', nargs='?', 
	help='output prefix for unalign and align mass fraction trajectory')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.2')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np

# default for args
args.omassf = args.output + '.mass'

## check vaulable setting
print("="*30)
if args.axis < 0 or args.axis > 2:
	raise ValueError("wrong input of axis for histogram")

## timer
start_proc, start_prof = hjung.time.init()

## read a topology and a trajectory using module MDAnalysis with selection
coordinates, unit_cells = hjung.io.read_trr_3d_select1(args.structure, args.input, args.select, 'pos')

## reduce 3d-coordinates to 1d-coordinates
unit_cells = hjung.array.convert_unitcell_3d(unit_cells, args.structure, args.input)
unit_cells_1d   = unit_cells[:,args.axis]
coordinates_1d = coordinates[:,:,args.axis]

## number histograms for each frame 
number_1d_t, bin_1d_t = hjung.analyze.histo_t_1d_nbin(coordinates_1d, unit_cells_1d, args.nbin) 
print("Done: making number trajectory with respect to bins")

## save number histogram trajectory
np.savetxt(args.omass, number_1d_t, 
	header='[%d, %d], mass1 fraction by molecules in nbins, %d' \
	%(len(number_1d_t),args.nbin,args.nbin), fmt='%f', comments='# ')
np.save(args.omass, number_1d_t)

## bin-size info 
box_axis_avg, box_axis_std = hjung.coord.box_1d_mode(unit_cells_1d,'box-z length','v')
print(" bin size = {0:.5f} +- {1:.5f}".format((box_axis_avg/float(args.nbin)),(box_axis_std/float(args.nbin))))

## timer
hjung.time.end_print(start_proc, start_prof)