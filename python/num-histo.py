#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 11/15/2016

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Number Probability of identically selected atoms in a trajectory for HJUNG module')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr structure file')
parser.add_argument('-select', '--select', nargs='?', 
	help='a file with a command-line for select_atoms in MDAnalysis')
parser.add_argument('-bin', '--bin', default=2.0, nargs='?', 
	help='bin size (A), defult 2A')
parser.add_argument('-axis', '--axis', default=2, nargs='?', 
	help='which axis for histogram (x axis (0), y axis (1), z axis (2))')
parser.add_argument('--nojump', action='store_true',
	help='re-positioning (x,y,z) within an unit cell (option)')
parser.add_argument('-onum', '--outputnumber', default='traj.nums', nargs='?', 
	help='output file for number trajectory')
parser.add_argument('-obin', '--outputbin', default='traj.bins', nargs='?', 
	help='output file for bin trajectory')
parser.add_argument('-o', '--output', default='traj.nprob', nargs='?', 
	help='output file for number probabilty')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# default for args
args.input = args.input if args.input is not None else 'traj.trr'
args.structure = args.structure if args.structure is not None else 'topol.tpr'
args.output = args.output if args.output is not None else 'traj.nprob'
args.outputnumber = args.outputnumber if args.outputnumber is not None else 'traj.nums'
args.outputbin = args.outputbin if args.outputbin is not None else 'traj.bins'
args.bin = args.bin if args.bin is not None else 2.0
args.axis = args.axis if args.axis is not None else 2
args.bin = float(args.bin)
args.axis = int(args.axis)

## Check arguments for log
print("===============================")
print("input filename   = ", args.input)
print("str filename     = ", args.structure)
if args.input is not None:
	print("select filename  = ", args.select)
print("No Jump if yes   = ", args.nojump)
print("bin size (A)     = ", args.bin)
print("axis [0:2]       = ", args.axis)
print("output number trajectory filename = ", args.outputnumber)
print("output bin trajectory filename = ", args.outputbin)
print("output number probability filename  = ", args.output)

## check vaulable setting
if args.bin < 1.0:
	raise ValueError("Too small bin size is usually not valid for histogram")
if args.axis < 0 or args.axis > 2:
	raise ValueError("wrong input of axis for histogram")

## timer
import time
start_clock = time.clock() # process time
start_wall = time.time() # wall time

## import hjung
import hjung
from hjung import *

## read a topology and a trajectory using module MDAnalysis with selection
print("===============================")
coordinates, unit_cells = hjung.io.read_coord_trr_3d(args.structure, args.input, args.select)
print("Done: reading trajectory and topology file")

## wrap positions within unit cell
if args.nojump:
	print("===============================")
	coordinates = hjung.coord.pbc_nojump_t(coordinates, unit_cells)
	print("Done: nojumping coordinates")

## reduce 3d-coordinates to 1d-coordinates
coordinates_1d = coordinates[:,:,args.axis]
unit_cells_1d = unit_cells[:,args.axis]

## number histograms for each frame 
print("===============================")
number_t_1d, bin_t_1d = hjung.analyze.histo_t_1d(coordinates_1d, unit_cells_1d, args.bin) 
print("Done: making number trajectory with respect to bins")

## save number histograms
print("===============================")
import numpy as np
np.savetxt(args.outputnumber, number_t_1d, fmt='%d')
np.savetxt(args.outputbin, bin_t_1d)
print("Finished saving files of number and bin trajectory")

## make an array for probability of numbers and histogram, then save the probability in a file
number_1d = hjung.array.merge_t_to_1d(number_t_1d)

## make 1d density probability of number
print("===============================")
prob_num_1d, prob_num_bin_1d = hjung.analyze.histo_num_dens_1d(number_1d)
print("Done: making a time- and bin-averaged probability of numbers")

## re-assign bin value for writing a file
print("===============================")
for i in range(1,len(prob_num_bin_1d)-1):
	prob_num_bin_1d[i-1] = 0.50*(prob_num_bin_1d[i-1]+prob_num_bin_1d[i])
prob_num_bin_1d = prob_num_bin_1d[0:len(prob_num_bin_1d)-1]
## save file
np.savetxt(args.output, np.column_stack((prob_num_bin_1d, prob_num_1d)), fmt="%7.3f %0.8f")
print("Finished saving files for probability of number")

## translate to get the lowest RMSE

## save the density profile trajectory

## block average

## save the block average

print("===============================")
print(time.clock() - start_clock, "seconds process time")
print(time.time() - start_wall, "seconds wall time")