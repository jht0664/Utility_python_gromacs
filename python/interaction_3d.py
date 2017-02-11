#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 1/03/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='vector trajectory of interactions of two selected atoms within cutoff distance')
## args
parser.add_argument('-b', '--b', nargs='?', help='begenning index of frame to analyze')
parser.add_argument('-e', '--e', nargs='?', help='end index of frame to analyze')
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr structure file')
parser.add_argument('-select1', '--select1', nargs='?', 
	help='a file with a command-line for first select_atoms in MDAnalysis')
parser.add_argument('-select2', '--select2', nargs='?', 
	help='a file with a command-line for second select_atoms in MDAnalysis')
parser.add_argument('-cutoff', '--cutoff', nargs='?', 
	help='cut-off distance (A), which should not exceed half of the smallest box length')
parser.add_argument('-plot3d', '--plot3d', default='traj.plot3d', nargs='?', 
	help='output file for 3d plot')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# default for args
args.b = args.b if args.b is not None else 0
args.b = int(args.b)
args.e = int(args.e)
args.input = args.input if args.input is not None else 'traj.trr'
args.structure = args.structure if args.structure is not None else 'topol.tpr'
args.plot3d = args.plot3d if args.plot3d is not None else 'traj.plot3d'
args.cutoff = float(args.cutoff)

## Check arguments for log
print("===============================")
print("begining iframe  = ", args.b)
print("end      iframe  = ", args.e)
print("input filename   = ", args.input)
print("str filename     = ", args.structure)
print("first select filename  = ", args.select1)
print("second select filename = ", args.select2)
print("cutoff (A)       = ", args.cutoff)
print("output 3d plot filename = ", args.plot3d)

## timer
import time
start_clock = time.clock() # process time
start_wall = time.time() # wall time

## import modules
import hjung
from hjung import *
import numpy as np

## read a topology and a trajectory using module MDAnalysis with selection
print("===============================")
group1_coord, box1 = hjung.io.read_coord_trr_3d(args.structure, args.input, args.select1)
group2_coord, box2 = hjung.io.read_coord_trr_3d(args.structure, args.input, args.select2)
print("Done: reading trajectory and topology file using select atoms")
if args.e is None:
	args.e = int(len(box1))

## check if the box is supported; rectangular
## assume unit_cells1 is the same as unit_cells2
if hjung.coord.rect_unit_cell(box1) is False:
	raise ValueError("First unit cell is not retangular which is not supported here")
## check if cutoff distance does not beyond the half of the smallest box length
if hjung.coord.dist_beyond_box(box1,args.cutoff) is True:
	raise ValueError("cutoff distance does not beyond the half of the smallest length of the first box")

## figure out interactions within the cut-off distance 
##  and the z-potion (center of the interaction distance) for each frame
##  assume # frames of unit_cells1 and unit_cells2 are the same
##  make and save number histogram for each frame
with open(args.plot3d, 'w') as output_file:
	output_file.write("## start points of vectors, direction of vector, magnitude(1/distance) \n")
	output_file.write("## x, y, z, dx, dy, dz, 1/r (unit, A) \n")

for iframe in range(args.b,args.e):
	print("current time step: %d" % iframe)
	result = hjung.analyze.rdf_cut_dt(group1_coord[iframe],group2_coord[iframe],box1[iframe],args.cutoff)
	file_append = open(args.plot3d,'ab')
	np.savetxt(file_append,result, header='time step: %d' % iframe,fmt='%g',comments='# ')
	file_append.close()
	# grouping dataset for gnuplot
	file_append = open(args.plot3d,'a')
	file_append.write('\n\n')
	file_append.close()

print("===============================")
print(time.clock() - start_clock, "seconds process time")
print(time.time() - start_wall, "seconds wall time")