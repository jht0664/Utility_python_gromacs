#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on  2/11/2017
# ver 0.2 - support pdb and dcd files for openmm on 5/8/2017
# ver 0.3 - support xtc trajectory files for Monte Carlo using "reduce_unitcells_3d_to_1d" on 6/6/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='2D Delta_number autocorrelation for selected molecules')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr or .gro structure file')
parser.add_argument('-select', '--select', nargs='?', 
	help='a file with a command-line for select_atoms in MDAnalysis')
parser.add_argument('-nbin', '--nbin', nargs='?', 
	help='number of bins')
parser.add_argument('-tol', '--tol', default=0.0, nargs='?', 
	help='tolerance for block average (> 0 and < 1). (recommend 1.0 for 1st trial). If 0, no block average. If > 1, # frames to average')
parser.add_argument('-axis', '--axis', default=2, nargs='?', 
	help='which axis for histogram (x axis (0), y axis (1), z axis (2))')
parser.add_argument('-o', '--output', default='traj.dnums', nargs='?', 
	help='output prefix for delta_number trajectory, spatial annd temporal delta_number autocorrelation, and delta_number trajectory after alignment using convolution')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.2')
# read args
args = parser.parse_args()
# default for args
args.input = args.input if args.input is not None else 'traj.trr'
args.structure = args.structure if args.structure is not None else 'topol.tpr'
args.odnum = args.output
args.odnumas = args.output + '.acfs'
args.odnumat = args.output + '.acft'
args.oalign = args.output + '.align'
args.axis = args.axis if args.axis is not None else 2
args.axis = int(args.axis)
args.nbin = int(args.nbin)
args.tol = args.tol if args.tol is not None else 0.0
args.tol = float(args.tol)


## Check arguments for log
print("===============================")
print("input filename   = ", args.input)
print("str filename     = ", args.structure)
if args.select is not None:
	print("select filename  = ", args.select)
print("number of bins   = ", args.nbin)
if args.tol == 0.0:
	print("Set no bloack average")
elif args.tol <= 1.0: 
	print("tolerance for block average = %f" %args.tol)
elif args.tol > 1.0:
	print("set block length = %d" %(int(args.tol)))
print("axis [0:2]       = ", args.axis)
print("output delta_number trajectory = ", args.odnum)
print("output spatial delta_number autocorrelation = ", args.odnumas)
print("output temporal delta_number autocorrelation = ", args.odnumat)
print("output aligned delta_number trajectory = ", args.oalign)

## check vaulable setting
if args.nbin < 10:
	raise ValueError("Too small nbin is usually not valid. (recommend > 10)")
if args.axis < 0 or args.axis > 2:
	raise ValueError("wrong input of axis for histogram")
if args.tol < 0.0:
	raise ValueError("wrong input of tolerance, %f" %args.tol)
elif args.tol >= 1.0: 
	print("Warning: tolerance %f is assigned to block_size, %d" %(args.tol, int(args.tol)))
else:
	print("="*30)

## timer
import time
start_clock = time.clock() # process time
start_wall = time.time() # wall time

## import modules
import hjung
from hjung import *
import numpy as np

## read a topology and a trajectory using module MDAnalysis with selection
print("="*30)
coordinates, unit_cells = hjung.io.read_coord_trr_3d(args.structure, args.input, args.select)
print("Done: reading trajectory and topology file")

## reduce 3d-coordinates to 1d-coordinates
coordinates_1d = coordinates[:,:,args.axis]
unit_cells_1d = hjung.array.reduce_unitcells_3d_to_1d(unit_cells, args.axis, args.structure, args.input)

## number histograms for each frame 
print("="*30)
number_t_1d, bin_t_1d = hjung.analyze.histo_t_1d_nbin(coordinates_1d, unit_cells_1d, args.nbin) 
print("Done: making number trajectory with respect to bins")

## block average to get stable volume (or number density)
block_length = 1
box_axis_avg, box_axis_std = hjung.coord.box_1d(unit_cells_1d)
print("box length avg = %f" %box_axis_avg)
print("bin size avg = %f" %(box_axis_avg/float(args.nbin)))
print("box length std = %f" %box_axis_std)
if args.tol > 0.0 and args.tol <= 1.0:
	print("="*30)
	print("To optimize, we use unit cell length on the axis you select.")
	block_length = hjung.analyze.opt_block_length_1d(unit_cells_1d,args.tol) 
	print("Done: optimize block length")
elif args.tol > 1: 
	block_length = int(args.tol)
if block_length > 1:
	print("="*30)
	number_t_1d = hjung.analyze.block_average_1d(number_t_1d,block_length) 

## save number histogram trajectory
#print("===============================")
#import numpy as np
#np.savetxt('traj.numb', number_t_1d, fmt='%f')
#print("Finished saving files of number and bin trajectory")

## Calculate delta_number, N[j-th frame, i-th slab] - <N>
# Assume: total #particle in a frame is fixed. (No particle incertion or deletion during simulation)
print("="*30)
ref_avg = np.mean(number_t_1d[0])
print("Assume: No particle incertion or deletion during simulation.")
delta_number_t_1d = number_t_1d - ref_avg
np.savetxt(args.odnum, delta_number_t_1d, 
	header='avg.number = %f, block_length = %d, generated by Hyuntae python code' %(ref_avg,block_length), fmt='%f', comments='# ')
print("Finished saving a file of delta_number trajectory.")

## autocorrelation function with periodic box system, but statistical autocorrelation on time trajectory
# \Delta N(i,{ t }_{ 1 })\ast \Delta N(j,{ t }_{ 2 })
# i.e. \Delta N(r)=N(r)-\left< c \right>  This is what we did previous step
# assume we use periodic condition for spatial coordinate.
print("="*30)
# set matrix spatial delta_delta_number 
print("spatial delta_number fluctuation")
acf_1d_wrap = hjung.analyze.autocorr_1d_t(delta_number_t_1d, 'wrap') 
slab_shift = int(len(acf_1d_wrap[0])/2.0)
np.savetxt(args.odnumas, acf_1d_wrap, 
	header='spatial autocorr(slab_lag,i_frame) for delta_number, Plot u ($1-%d):2:3 when block_length = %d'	 
	%(slab_shift,block_length), fmt='%f', comments='# ')
# set matrix temporal delta_delta_number 
#print("temporal delta_number fluctuation")
#transp_delta_number = delta_number_t_1d.transpose()
#acf_1d_temp = hjung.analyze.autocorr_1d_t(transp_delta_number, 'constant') 
#time_shift = int(len(acf_1d_temp[0])/2.0)
#np.savetxt(args.odnumat, acf_1d_temp, 
#	header='Temporal autocorr(time_lag,i_slab) for delta_number. Plot u ($1-%d):2:3 when block_length = %d' 
#	%(time_shift,block_length), fmt='%f', comments='# ')
#print("Finished saving a file of spatial and temporal delta_number autocorrelation")

# make gnuplot script
print("="*30)
print("for gnuplotting, please use following code:")
print("load '/home/hjung52/Utility/gnuplot/palette.plt'")
print("set pm3d map")
print("set multiplot layout 2, 2 title 'PEO plotting'")
#for first plot
print("set tmargin 2")
print("unset key")
print("set title 'PEO delta_Number'")
print("splot 'peo.dnum' u ($1*bin_length):($2*traj_1frame_ns*%d):($3/30) matrix" % (block_length))
# for second plot
print("set title 'PEO spatial autocorrelation with PBC'")
# load '/home/hjung52/Utility/gnuplot/palette.plt'
# set pm3d map
# set multiplot layout 2, 2 title 'PEO plotting'
# set tmargin 2
# unset key
# set title 'PEO delta_Number'
# set cbr[-2.5:2.5]
# splot 'peo.dnum' u ($1*0.19888):($2*0.01*2):($3/30) matrix
# set title 'PEO spatial autocorrelation with PBC'
# set cbr[-1:1]
# splot 'peo.dnumas' u ($1*0.19888-12.9274/2+0.09888):($2*0.02):3 matrix
# splot 'peo.dnumat' u ($1*0.02-10):($2*0.19888):3 matrix
# unset multiplot

## align delta_number using spatial autocorrelation function
print("="*30)
align_delta_number_t_1d =  hjung.analyze.align_acf(delta_number_t_1d, acf_1d_wrap, 'wrap') 
np.savetxt(args.oalign, align_delta_number_t_1d, 
	header='aligned delta_number by hjung.analyze.align_acf function', fmt='%f', comments='# ')
print("Finished saving a file of aligned delta_number using acf")


print("="*30)
print(time.clock() - start_clock, "seconds process time")
print(time.time() - start_wall, "seconds wall time")