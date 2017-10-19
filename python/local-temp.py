#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 10/18/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculation of local temperatures of binary mixture in slab geometry')
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
parser.add_argument('-o', '--output', default='traj', nargs='?', 
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

## default for args output
# 3d xyz
args.otmpd = args.output + '.temp'            # based temperature from delta_dist
args.otmpda = args.output + '.temp.algin'     # corrected temperature
args.otmpv = args.output + '.temp.sys'       # based temperature from velocity trajectory
args.otmpva = args.output + '.temp.sys.align' # corrected temperature?
# 1d z
args.otmpdz = args.output + '.temp.z'            # based temperature from delta_dist
args.otmpdaz = args.output + '.temp.z.algin'     # corrected temperature
args.otmpvz = args.output + '.temp.z.sys'       # based temperature from velocity trajectory
args.otmpvaz = args.output + '.temp.z.sys.align' # corrected temperature?
# 2d xy
args.otmpdxy = args.output + '.temp.xy'            # based temperature from delta_dist
args.otmpdaxy = args.output + '.temp.xy.algin'     # corrected temperature
args.otmpvxy = args.output + '.temp.xy.sys'       # based temperature from velocity trajectory
args.otmpvaxy = args.output + '.temp.xy.sys.align' # corrected temperature?

## check arguments for log
print("===============================")
print("input filename   = ", args.input)
print("str filename     = ", args.structure)
print("mass info filename = ", args.mass)
print("select1 filename  = ", args.select1)
print("select2 filename  = ", args.select2)
hjung.blockavg.print_init(args.tol)
print("number of bins   = ", args.nbin)
print("axis [0:2]       = ", args.axis)
print("output temp. filenames = ", args.otmpd, args.otmpda, args.otmpv, args.otmpva)

## check vaulable setting
print("="*30)
if args.axis < 0 or args.axis > 2:
	raise ValueError("wrong input of axis for histogram")
args.tol = hjung.blockavg.check(args.tol)
print("="*30)

## timer
start_proc, start_prof = hjung.time.init()

## read a topology and a trajectory using module MDAnalysis with selection
print("="*30)
data1, data2, unit_cells = hjung.io.read_trr_3d_select2(args.structure, args.input, args.select1, args.select2, 'pos vel')
coordinates1 = data1[0] # unit = A
velocity1 = data1[1] # unit = A/ps
coordinates2 = data2[0] 
velocity2 = data2[1]
print("Done: reading trajectory and topology file")
#for ivel in velocity1: 
#	for iatom in ivel:
#		np.square(iatom)

## reduce 3d-coordinates to 1d-coordinates
# complicate way?
#coordinates1_1d = hjung.array.reduce_3d_to_1d(coordinates1,args.axis)
#velocity1_1d    = hjung.array.reduce_3d_to_1d(velocity1,args.axis)
#coordinates2_1d = hjung.array.reduce_3d_to_1d(coordinates2,args.axis)
#velocity2_1d    = hjung.array.reduce_3d_to_1d(velocity2,args.axis)
coordinates1_1d = coordinates1[:,:,args.axis]
velocity1_1d = velocity1[:,:,args.axis]
coordinates2_1d = coordinates2[:,:,args.axis]
velocity1_1d = velocity1[:,:,args.axis]
unit_cells_1d = hjung.array.reduce_unitcells_3d_to_1d(unit_cells, args.axis, args.structure, args.input)

# unit cell info
box_axis_avg, box_axis_std = hjung.coord.box_1d(unit_cells_1d)
print(" box length avg = %f" %box_axis_avg)
print(" bin size avg = %f" %(box_axis_avg/float(args.nbin)))
print(" box length std = %f" %box_axis_std)

## number histograms for each frame 
print("="*30)
number1_1d_t, bin_1d_t = hjung.analyze.histo_t_1d_nbin(coordinates1_1d, unit_cells_1d, args.nbin) 
number2_1d_t, bin_1d_t = hjung.analyze.histo_t_1d_nbin(coordinates2_1d, unit_cells_1d, args.nbin) 
n_frames = len(coordinates1_1d)
n_atoms = len(coordinates1_1d[0])
index_coord1 = np.zeros((n_frames,n_atoms),dtype=int)
index_coord2 = np.zeros((n_frames,n_atoms),dtype=int)
for i_frame in range(n_frames):
	index_coord1[i_frame] = np.digitize(coordinates1_1d[i_frame],bin_1d_t[i_frame], right=False)
	index_coord2[i_frame] = np.digitize(coordinates2_1d[i_frame],bin_1d_t[i_frame], right=False)
# shift from ibin[1:nbin+1] to bin index[0:nbin]
index_coord1 = index_coord1 - 1
index_coord2 = index_coord2 - 1
# if bin index = nbin, it should belong 0-th, not nbin bin by pbc
index_coord1 = np.mod(index_coord1,args.nbin)
index_coord2 = np.mod(index_coord2,args.nbin)
print("Done: making number trajectory with respect to bins")

## 3d temperature
print("="*30)
print("System temperature calculation")
n_frames = len(index_coord1)
temp_sys_3d = np.zeros((n_frames,args.nbin))
temp_sys_1d = np.zeros((n_frames,args.nbin))
temp1_sys_3d = np.zeros((n_frames,args.nbin))
temp1_sys_1d = np.zeros((n_frames,args.nbin))
temp2_sys_3d = np.zeros((n_frames,args.nbin))
temp2_sys_1d = np.zeros((n_frames,args.nbin))
mod_frame = hjung.io.process_init()
for i_frame in range(n_frames):
	i_coord1 = index_coord1[i_frame]
	i_coord2 = index_coord2[i_frame]
	n_atoms1 = len(i_coord1)
	count_natoms1_bin = np.zeros(args.nbin)
	count_natoms2_bin = np.zeros(args.nbin)
	# add v**2 data for select1
	for i_atom in range(n_atoms1):
		i_bin = int(i_coord1[i_atom])
		print(i_bin)
		count_natoms1_bin[i_bin] = count_natoms1_bin[i_bin] + 1.0
		vxyz = velocity1[i_frame][i_atom]
		temp1_sys_3d[i_frame][i_bin] = temp1_sys_3d[i_frame][i_bin] + np.sum(np.square(vxyz))/3.0 
		vz = velocity1[i_frame][i_atom][2]
		temp1_sys_1d[i_frame][i_bin] = temp1_sys_1d[i_frame][i_bin] + np.square(vz)/3.0
	# average v**2
	for i_bin in range(args.nbin):
		if count_natoms1_bin[i_bin] == 0:
			raise ValueError("too small bin length. Decrease nbin!")
		temp1_sys_3d[i_frame][i_bin] = temp1_sys_3d[i_frame][i_bin]/count_natoms1_bin[i_bin]
		temp1_sys_1d[i_frame][i_bin] = temp1_sys_1d[i_frame][i_bin]/count_natoms1_bin[i_bin]

	# add v**2 data for select2
	n_atoms2 = len(i_coord2)
	for i_atom in range(n_atoms2):
		i_bin = int(i_coord2[i_atom])
		count_natoms2_bin[i_bin] = count_natoms2_bin[i_bin] + 1.0
		vxyz = velocity2[i_frame][i_atom]
		temp2_sys_3d[i_frame][i_bin] = temp2_sys_3d[i_frame][i_bin] + np.sum(np.square(vxyz))/3.0
		vz = velocity2[i_frame][i_atom][2]
		temp2_sys_1d[i_frame][i_bin] = temp2_sys_1d[i_frame][i_bin] + np.square(vz)/3.0
	# average v**2
	for i_bin in range(args.nbin):
		if count_natoms2_bin[i_bin] == 0:
			raise ValueError("too small bin length. Decrease nbin!")
		temp2_sys_3d[i_frame][i_bin] = temp2_sys_3d[i_frame][i_bin]/count_natoms2_bin[i_bin]
		temp2_sys_1d[i_frame][i_bin] = temp2_sys_1d[i_frame][i_bin]/count_natoms2_bin[i_bin]
	# print process
	mod_frame = hjung.io.process_print(i_frame+1, n_frames, mod_frame)
temp_sys_3d = temp1_sys_3d + temp2_sys_3d
temp_sys_1d = temp1_sys_1d + temp2_sys_1d
temp_sys_2d = temp_sys_3d - temp_sys_1d
i_frame = n_frames - 1 
print("{} frame -> sum T_xyz {}, T_z {}, T_xy {}".format(i_frame,np.sum(temp_sys_3d[i_frame]),np.sum(temp_sys_1d[i_frame]),np.sum(temp_sys_2d[i_frame])))
i_frame = n_frames - 2
print("{} frame -> sum T_xyz {}, T_z {}, T_xy {}".format(i_frame,np.sum(temp_sys_3d[i_frame]),np.sum(temp_sys_1d[i_frame]),np.sum(temp_sys_2d[i_frame])))
print("Done: 3d temperature of system")
np.savetxt(args.otmpv, temp_sys_3d, 
	header='[%d, %d], total temperature (xyz) in nbins, %d' \
	%(len(temp_sys_3d),args.nbin,args.nbin), fmt='%f', comments='# ')
np.savetxt(args.otmpvxy, temp_sys_2d, 
	header='[%d, %d], temperature (xy) in nbins, %d' \
	%(len(temp_sys_2d),args.nbin,args.nbin), fmt='%f', comments='# ')
np.savetxt(args.otmpvz, temp_sys_1d, 
	header='[%d, %d], temperature (z) in nbins, %d' \
	%(len(temp_sys_1d),args.nbin,args.nbin), fmt='%f', comments='# ')
# difference from average
diff_sys_3d = np.zeros((n_frames,args.nbin))
diff_sys_2d = np.zeros((n_frames,args.nbin))
diff_sys_1d = np.zeros((n_frames,args.nbin))
for i_frame in range(n_frames):
	diff_sys_3d[i_frame] = temp_sys_3d[i_frame] - np.average(temp_sys_3d[i_frame])
	diff_sys_2d[i_frame] = temp_sys_2d[i_frame] - np.average(temp_sys_2d[i_frame])
	diff_sys_1d[i_frame] = temp_sys_1d[i_frame] - np.average(temp_sys_1d[i_frame])
np.savetxt('diff3', diff_sys_3d, 
	header='[%d, %d], diff total temperature (xyz) in nbins, %d' \
	%(len(diff_sys_3d),args.nbin,args.nbin), fmt='%f', comments='# ')
np.savetxt('diff2', diff_sys_2d, 
	header='[%d, %d], diff temperature (xy) in nbins, %d' \
	%(len(diff_sys_2d),args.nbin,args.nbin), fmt='%f', comments='# ')
np.savetxt('diff1', diff_sys_1d, 
	header='[%d, %d], diff temperature (z) in nbins, %d' \
	%(len(diff_sys_1d),args.nbin,args.nbin), fmt='%f', comments='# ')

# obtain alignment shift
print("="*30)
print("Find alignment shift")
## read mass info
mw, divider = hjung.io.read_mass2(args.mass)
## Calculate mass fraction of each bins
mass1_1d_t = number1_1d_t*mw[0]/divider[0]
mass2_1d_t = number2_1d_t*mw[1]/divider[1]
totalmass_1d_t = np.array(mass1_1d_t + mass2_1d_t,dtype=np.float)
massfrac_1d_t = np.divide(mass1_1d_t,totalmass_1d_t)
## Align mass fractions using autocorrelation function
acf_1d_t_wrap = hjung.analyze.autocorr_1d_t(massfrac_1d_t, 'wrap') 
align_shift = hjung.analyze.convolve_1d_t(acf_1d_t_wrap, massfrac_1d_t, 'wrap', 'max') 
print(" Convolution std = {}".format(np.std(align_shift)))
# shifting
for i_frame in range(len(temp_sys_3d)):
	shift_array3 = temp_sys_3d[i_frame]
	shift_array2 = temp_sys_2d[i_frame]
	shift_array1 = temp_sys_1d[i_frame]
	temp_sys_1d[i_frame] = np.roll(shift_array1, align_shift[i_frame]) 
	temp_sys_2d[i_frame] = np.roll(shift_array2, align_shift[i_frame])
	temp_sys_3d[i_frame] = np.roll(shift_array3, align_shift[i_frame])
print("Done: align 3d temperature of system")
np.savetxt(args.otmpva, temp_sys_3d, 
	header='[%d, %d], aligned total temperature (xyz) in nbins, %d' \
	%(len(temp_sys_3d),args.nbin,args.nbin), fmt='%f', comments='# ')
np.savetxt(args.otmpvaxy, temp_sys_2d, 
	header='[%d, %d], aligned temperature (xy) in nbins, %d' \
	%(len(temp_sys_2d),args.nbin,args.nbin), fmt='%f', comments='# ')
np.savetxt(args.otmpvaz, temp_sys_1d, 
	header='[%d, %d], aligned temperature (z) in nbins, %d' \
	%(len(temp_sys_1d),args.nbin,args.nbin), fmt='%f', comments='# ')



## timer
hjung.time.end_print(start_proc, start_prof)