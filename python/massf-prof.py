#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 02/25/2017
# ver 0.2 - support pdb and dcd files for openmm on 5/8/2017
# ver 0.3 - support xtc trajectory files for Monte Carlo using "reduce_unitcells_3d_to_1d" on 6/6/2017
# ver 0.4 - support block average module on 6/26/2017
# ver 0.5 - add some function: printing arguments and reduce lines as for default setting

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='1D mass fraction profile using counting number of particles')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr structure file')
parser.add_argument('-m', '--mass', nargs='?', 
	help='divider for normalization and masses for selected molecules')
parser.add_argument('-select1', '--select1', nargs='?', 
	help='a file1 with a command-line for select_atoms in MDAnalysis')
parser.add_argument('-select2', '--select2', nargs='?', 
	help='a file2 with a command-line for select_atoms in MDAnalysis')
parser.add_argument('-nbin', '--nbin', nargs='?', type=int,
	help='number of bins')
parser.add_argument('-tol', '--tol', default=0.0, nargs='?', type=float,
	help='tolerance for block average (> 0 and < 1). (recommend 1.0 for 1st trial). If 0, no block average. If > 1, # frames to average')
parser.add_argument('-axis', '--axis', default=2, nargs='?', type=int,
	help='which axis for histogram (x axis (0), y axis (1), z axis (2))')
parser.add_argument('-align', '--align', default='YES', nargs='?',
	help='Run alignment? or not? (YES/NO)')
parser.add_argument('-o', '--output', default='traj.massf', nargs='?', 
	help='output prefix for unalign and align mass1 fraction trajectory')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.4')
# read args
args = parser.parse_args()
# check args
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np

# default for args
args.oalign = args.output + '.align'

## Check arguments for log
print("===============================")
print("input filename   = ", args.input)
print("str filename     = ", args.structure)
print("mass info filename = ", args.mass)
print("select1 filename  = ", args.select1)
print("select2 filename  = ", args.select2)
hjung.blockavg.print_init(args.tol)
print("number of bins   = ", args.nbin)
print("axis [0:2]       = ", args.axis)
print("output mass frac filename  = ", args.output)
print("output align mass frac filename  = ", args.oalign)

## check vaulable setting
if args.axis < 0 or args.axis > 2:
	raise ValueError("wrong input of axis for histogram")
args.tol = hjung.blockavg.check(args.tol)

## timer
start_proc, start_prof = hjung.time.init()

## read a topology and a trajectory using module MDAnalysis with selection
print("="*30)
coordinates1, coordinates2, unit_cells = hjung.io.read_coord_trr_3d_select2(args.structure, args.input, args.select1, args.select2)
print("Done: reading trajectory and topology file")

## reduce 3d-coordinates to 1d-coordinates
coordinates1_1d = coordinates1[:,:,args.axis]
coordinates2_1d = coordinates2[:,:,args.axis]
unit_cells_1d = hjung.array.reduce_unitcells_3d_to_1d(unit_cells, args.axis, args.structure, args.input)

box_axis_avg, box_axis_std = hjung.coord.box_1d(unit_cells_1d)
print(" box length avg = %f" %box_axis_avg)
print(" bin size avg = %f" %(box_axis_avg/float(args.nbin)))
print(" box length std = %f" %box_axis_std)

## number histograms for each frame 
print("="*30)
number1_1d_t, bin_1d_t = hjung.analyze.histo_t_1d_nbin(coordinates1_1d, unit_cells_1d, args.nbin) 
number2_1d_t, bin_1d_t = hjung.analyze.histo_t_1d_nbin(coordinates2_1d, unit_cells_1d, args.nbin) 
print("Done: making number trajectory with respect to bins")

## read args.mass file
print("="*30)
try:
	massinfo = open(args.mass, 'r')
except IOError:
	print("Problem with opening ",args.mass)
	exit()
divider = []
mw = []
for line in massinfo:
	line = line.strip()
	line_m = line.rsplit()
	divider.append(float(line_m[0]))
	mw.append(float(line_m[1]))
massinfo.close()
divider = np.array(divider)
mw = np.array(mw)
if len(divider) != 2 or len(mw) != 2:
	ValueError("Wrong format in %s file" %args.mass)
print("dividers[select1,select2] = %s" %divider)
print("mw[select1,select2] = %s" %mw)
## Calculate mass fraction of each bins
mass1_1d_t = number1_1d_t*mw[0]/divider[0]
#print("mass1_1d_t %s" %mass1_1d_t)
mass2_1d_t = number2_1d_t*mw[1]/divider[1]
#print("mass2_1d_t %s" %mass2_1d_t)
massfrac_1d_t = mass1_1d_t/(mass1_1d_t+mass2_1d_t)
#print("massfrac_1d_t %s" %massfrac_1d_t)

## block average
print("="*30)
massfrac_1d_t, block_length = hjung.blockavg.main_1d(massfrac_1d_t, unit_cells_1d, args.tol) 

## save number histogram trajectory
np.savetxt(args.output, massfrac_1d_t, 
	header='[%d, %d], mass1 fraction by molecules in nbins, %d' \
	%(len(mass1_1d_t),args.nbin,args.nbin), fmt='%f', comments='# ')

## Align mass fractions using autocorrelation function
acf_1d_t_wrap = hjung.analyze.autocorr_1d_t(massfrac_1d_t, 'wrap') 
if (args.align == 'YES'):
	align_massfrac_1d_t =  hjung.analyze.align_acf(massfrac_1d_t, acf_1d_t_wrap, 'wrap') 
else:
	align_massfrac_1d_t = massfrac_1d_t
np.savetxt(args.oalign, align_massfrac_1d_t, 
	header='%d, %d, aligned mass1 fraction by ACF and molecules in nbins' \
	%(len(mass1_1d_t),args.nbin), fmt='%f', comments='# ')
print("Finished saving unalign and align output files")

## timer
hjung.time.end_print(start_proc, start_prof)