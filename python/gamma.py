#!/usr/bin/env python3
# ver 0.1 - taking from massf-prof.py and coding python by Hyuntae Jung on 11/19/2017 

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='surface tension calculation using test volume method for WindoM-Rowlinson model')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='structure file')
parser.add_argument('-dv', '--dv', nargs='?', type=float,
	help='decrease of del_volume (should be absolute value)')
parser.add_argument('-nw', '--nwindow', nargs='?', type=int,
	help='number of windows for decreasing del_volume')
parser.add_argument('-select1', '--select1', nargs='?', 
	help='a file1 with a command-line for select_atoms in MDAnalysis')
parser.add_argument('-select2', '--select2', nargs='?', 
	help='a file2 with a command-line for select_atoms in MDAnalysis')
parser.add_argument('-o', '--output', default='traj', nargs='?', 
	help='output prefix for unalign and align mass1 fraction trajectory')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# check args
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np
import copy 

# default for args
args.op = args.output + '.pres'
args.opx = args.op + '.x'
args.opy = args.op + '.y'
args.opz = args.op + '.z'
args.os = args.output + '.surf'

## timer
start_proc, start_prof = hjung.time.init()

## read a topology and a trajectory using module MDAnalysis with selection
print("="*30)
coordinates1, coordinates2, unit_cells = hjung.io.read_trr_3d_select2(args.structure, args.input, args.select1, args.select2, 'pos')
print("Done: reading trajectory and topology file")

## convert unit_cells in universal format
unit_cells = hjung.array.convert_unitcell_3d(unit_cells, args.structure, args.input)

## calculate scale_length from del_volume
scale_length = np.zeros(args.nwindow)
del_volume = -args.dv 
for i in range(args.nwindow):
    dvol = 1.0 + del_volume*(i+1)
    if dvol < 0:
        raise ValueError("too big value of args.dv and/or many windows, leading to {}".foramt(dvol))
    scale_length[i] = dvol**(1.0/3.0)

# calculate number of overlapping
print("="*30)
n_frames = len(coordinates1)
n_atoms1 = len(coordinates1[0])
n_atoms2 = len(coordinates2[0])
n_overlap = np.zeros((n_frames,args.nwindow,3),dtype=int)
cell_xyz = unit_cells[0] # we assume NVT ensembles
mod_frame = hjung.io.process_init() # initialize print process
for i_frame in range(n_frames):
	for i in range(n_atoms1):
		print(i)
		for j in range(n_atoms2):
			dx = np.zeros(3)
			dx1 = np.zeros(3)
			for iaxis in range(3):
				dx[iaxis] = coordinates1[i_frame][i][iaxis] - coordinates1[i_frame][j][iaxis]
				dx1[iaxis] = dx[iaxis]/cell_xyz[iaxis]
				dx[iaxis] = dx[iaxis] - cell_xyz[iaxis]*np.around(dx1[iaxis])
			for iw in range(args.nwindow):
				for iaxis in range(3):
					dx_iaxis = copy.copy(dx)
					dx_iaxis[iaxis] = dx_iaxis[iaxis]*scale_length[iw]
					dist = np.sum(np.square(dx_iaxis))
					if dist < 1.0:
						n_overlap[i_frame][iw][iaxis] = n_overlap[i_frame][iw][iaxis] + 1
	mod_frame = hjung.io.process_print(i_frame+1, n_frames, mod_frame)	# print process
print("Done: calculate number of overlapping")

# calculate pressure tensor
print("="*30)
vol = cell_xyz[0]*cell_xyz[1]*cell_xyz[2]
t_atoms = n_atoms1 + n_atoms2
output_pres = np.zeros((n_frames,args.nwindow,3))
mod_frame = hjung.io.process_init() # initialize print process
for i_frame in range(n_frames):
	for iw in range(args.nwindow):
		for iaxis in range(3):
			output_pres[i_frame][iw][iaxis] = t_atoms/vol 
			+ n_overlap[i_frame][iw][iaxis]/((1.0-scale_length[iw])*vol)
	mod_frame = hjung.io.process_print(i_frame+1, n_frames, mod_frame)	# print process
output_surf = np.zeros((n_frames,args.nwindow))
for i_frame in range(n_frames):
	for iw in range(args.nwindow):
		output_surf[i_frame][iw] = cell_xyz[2]*(output_pres[i_frame][iw][2]-(output_pres[i_frame][iw][0]+output_pres[i_frame][iw][1])/2.0)
print("Done: calculate pressure tensor and surface tension")

## output
np.savetxt(args.os, output_surf, 
	header='[%d, %d], surface tension, %d' \
	%(n_frames,args.nwindown), fmt='%f', comments='# ')
np.savetxt(args.opx, output_pres[:][:][0], 
	header='[%d, %d], pressure_xx tensor, %d' \
	%(n_frames,args.nwindown), fmt='%f', comments='# ')
np.savetxt(args.opy, output_pres[:][:][1], 
	header='[%d, %d], pressure_yy tensor, %d' \
	%(n_frames,args.nwindown), fmt='%f', comments='# ')
np.savetxt(args.opz, output_pres[:][:][2], 
	header='[%d, %d], pressure_zz tensor, %d' \
	%(n_frames,args.nwindown), fmt='%f', comments='# ')

## timer
hjung.time.end_print(start_proc, start_prof)