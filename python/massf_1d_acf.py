#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 04/06/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='autocorrelation functions along time for a given position')
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
parser.add_argument('-nxbin', '--nxbin', default=-1, nargs='?', type=int,
	help='number of bins along x,y axis (or others except args.axis). If negative, use a length close to args.nzbin')
parser.add_argument('-nzbin', '--nzbin', nargs='?', type=int,
	help='number of bins along z axis (or args.axis) which should be the same when you got a.massf.opt.npy in massf_1d.py')
parser.add_argument('-axis', '--axis', default=2, nargs='?', type=int,
	help='which axis for histogram (x axis (0), y axis (1), z axis (2))')
parser.add_argument('-irem', '--irem', default='a.massf.remove.npy', nargs='?', 
	help='index of frames to remove from conv_1d.py')
parser.add_argument('-o', '--output', default='traj', nargs='?', 
	help='output prefix for unalign and 3D align mass1 fraction trajectory')
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
import numpy as np
from scipy import ndimage

# default for args
args.omassf = args.output + '.int_massf'
args.otmass = args.output + '.int_tmass'

## check vaulable setting
if args.axis < 0 or args.axis > 2:
	raise ValueError("wrong input of axis for histogram")
if args.axis != 2:
	raise ValueError("if the access")
## timer
start_proc, start_prof = hjung.time.init()

## read a topology and a trajectory using module MDAnalysis with selection
coordinates1, coordinates2, unit_cells = hjung.io.read_trr_3d_select2(args.structure, args.input, args.select1, args.select2, 'pos')
unit_cells = hjung.array.convert_unitcell_3d(unit_cells, args.structure, args.input)

## load other files to remove and shift
if 'NONE' in args.irem:
	print("you select no removing frames")
else:
	multilayer_iframes = np.load(args.irem)
	## remove time frames
	coordinates1 = np.delete(coordinates1,multilayer_iframes,axis=0)
	coordinates2 = np.delete(coordinates2,multilayer_iframes,axis=0)
	unit_cells = np.delete(unit_cells,multilayer_iframes,axis=0)

## determine nxbin from nzbin
if args.nxbin < 0:
	nxbin_length = np.mean(unit_cells[:,args.axis])/float(args.nzbin)
	unit_cells_x = np.mean(unit_cells[:,args.axis-1])
	args.nxbin = int(unit_cells_x/nxbin_length) # assume x,y dimension of unit_cells are same
binsize = np.array([args.nxbin,args.nxbin,args.nzbin])
nxbin_length = unit_cells[:,args.axis-1]/float(args.nxbin)
nzbin_length = unit_cells[:,args.axis]/float(args.nzbin)
print("we set nxbin = {} of which length (avg) = {}".format(args.nxbin,np.mean(nxbin_length)))
print("we set nzbin = {} of which length (avg) = {}".format(args.nzbin,np.mean(nzbin_length)))
nbin_length = np.column_stack((nxbin_length,nxbin_length,nzbin_length))

## generate coordinate with repect to bin_3d
nbin_length_reshape = nbin_length.reshape(-1,1,3)
coordinates1_bin_3d = (np.divide(coordinates1,nbin_length_reshape)).astype(int)
coordinates2_bin_3d = (np.divide(coordinates2,nbin_length_reshape)).astype(int)
#check wrong values
if np.amin(coordinates1_bin_3d) < 0:
	raise RuntimeError("negative values in axis {}?".format(np.amin(coordinates1_bin_3d)))
if np.max(np.amax(np.amax(coordinates1_bin_3d,axis=1),axis=0) > binsize): 
	raise RuntimeError("select1 beyond pbc like {} in ref {}?".format(np.amax(np.amax(coordinates1_bin_3d,axis=1),axis=0), binsize))
if np.max(np.amax(np.amax(coordinates2_bin_3d,axis=1),axis=0) > binsize): 
	raise RuntimeError("select2 beyond pbc like {} in ref {}?".format(np.amax(np.amax(coordinates2_bin_3d,axis=1),axis=0), binsize))

def histo_nbin(x, n_bins):
	# make histogram trajectory
	i_frame = 0
	if (np.amax(x) > n_bins) or (np.amin(x) < 0.0):
		print(" ##### {} {} {} Be aware of maximum value over (or less) box size, due to {} > {}, or {} < {} ##### ".format(i_frame,np.argmax(x),np.argmin(x),np.amax(x),n_bins,np.amin(x),0.0))
	return np.histogram(x, bins=n_bins, range=(0,n_bins))

def acf_1d_fn(data_1d, setmode):
	acf_data_1d = np.zeros(np.shape(data_1d))
	delta_data_1d = data_1d - data_1d.mean()
	acf_data_1d = ndimage.correlate(delta_data_1d,delta_data_1d,mode=setmode)
	acf_data_1d /= (data_1d.var()*len(delta_data_1d)) # normalize
	return acf_data_1d

## read args.mass file for weights
mw, divider = hjung.io.read_mass2(args.mass)

## make massf, do conv. and align
n_frames = len(coordinates1_bin_3d)
intrinsic_mf_profile = np.zeros((n_frames, args.nzbin))
intrinsic_tm_profile = np.zeros((n_frames, args.nzbin))
np.seterr(divide='ignore')
n_vacant_cells = 0
n_vacant_times = 0
n_vacant_avg = 0
n_vacant_avg_times = 0
TF_skip_frame = False
skip_frame_list = []
imod = hjung.time.process_init()
for iframe in range(n_frames):
	n_vacant_cells = 0
	n_vacant_times = 0				
	if TF_skip_frame:
		TF_skip_frame = False
	for ix in range(args.nxbin):
		if TF_skip_frame:
			break
		for iy in range(args.nxbin):
			if TF_skip_frame:
				break
			## number histograms for each frame 
			#print(" iframe = {}, ix = {}, iy = {}".format(iframe,ix,iy))
			#print(coordinates1_bin_3d)
			#print(coordinates1_bin_3d.shape)
			#print(np.where((coordinates1_bin_3d[iframe,:,0] == ix)))
			#print(coordinates1_bin_3d[iframe][np.where((coordinates1_bin_3d[iframe,:,0] == ix))])
			#print(np.where((coordinates1_bin_3d[iframe,:,0] == ix) & (coordinates1_bin_3d[iframe,:,1] == iy)))
			coord1_bin_1d_z = (coordinates1_bin_3d[iframe][np.where((coordinates1_bin_3d[iframe,:,0] == ix) & (coordinates1_bin_3d[iframe,:,1] == iy))])[:,2]
			coord2_bin_1d_z = (coordinates2_bin_3d[iframe][np.where((coordinates2_bin_3d[iframe,:,0] == ix) & (coordinates2_bin_3d[iframe,:,1] == iy))])[:,2]
			number1_1d, bin_1d = histo_nbin(coord1_bin_1d_z, args.nzbin) 
			number2_1d, bin_1d = histo_nbin(coord2_bin_1d_z, args.nzbin) 
			#print(number1_1d)
			#print(number2_1d)
			## Calculate mass fraction of each bins with weights
			mass1_1d = np.array(number1_1d*mw[0]/divider[0],dtype=np.float)
			mass2_1d = np.array(number2_1d*mw[1]/divider[1],dtype=np.float)
			totalmass_1d = mass1_1d + mass2_1d
			#print(int(np.sum(totalmass_1d)))
			#print(totalmass_1d)
			#print(number1_1d)
			#print(len(totalmass_1d))
			massfrac_1d = np.divide(mass1_1d,totalmass_1d)
			#print(massfrac_1d)
			#print(np.where(totalmass_1d == 0))
			list_vacant_cell = np.where(totalmass_1d == 0)[0] # or np.where(totalmass_1d == 0)[0])
			#print(list_vacant_cell)
			if len(list_vacant_cell) > 0:
				#print(np.diff(list_vacant_cell))
				#print(np.where(np.diff(list_vacant_cell) == 1)[0])
				#check if all neighbor cells has at least one particle to guess values in the vacancy cells
				if len(np.where(np.diff(list_vacant_cell) == 1)[0]) > 0:
					print("even neighbor cell does not have values and is vacant at {} frame.".format(iframe))
					TF_skip_frame = True
					skip_frame_list.append(iframe)
					#print(skip_frame_list)
					break
					#raise RuntimeError("even neighbor cell does not have values and is vacant.")
				#print("we estimate mean value to fill in {} vacant cells".format(len(list_vacant_cell)))
				n_vacant_cells = n_vacant_cells + len(list_vacant_cell)
				n_vacant_times = n_vacant_times + 1
				left_neighbor = np.remainder(list_vacant_cell - 1,args.nzbin)
				right_neighbor = np.remainder(list_vacant_cell + 1,args.nzbin)
				massfrac_1d[list_vacant_cell] = 0.50*(massfrac_1d[left_neighbor] + massfrac_1d[right_neighbor])
				#print(massfrac_1d)
				#raise ValueError("stop")			
			## autocorrelation
			#acf_1d = acf_1d_fn(massfrac_1d,'wrap')
			## check multilayers
			#step_1d = np.where(acf_1d < 0.0, -1., 1.)
			##print(acf_1d)
			##print(step_1d)
			#step_diff = np.diff(step_1d) # difference of neighbor element
			#step_n_up = int((step_diff > 0.0).sum())
			#step_n_down = int((step_diff < 0.0).sum())
			## the sum of number of positives or negatives should have 1 if phase separation noramlly occurs
			#if (step_n_up != 1) or (step_n_down != 1): 
			#	#print(" multi-domain problem ({} {}) at {}".format(step_n_up,step_n_down,iframe))
			#	TF_skip_frame = True
			#	skip_frame_list.append(iframe)
			#	#print(skip_frame_list)
			#	break
			step_1d = np.where(np.roll(np.arange(args.nzbin), int(args.nzbin/4)) < args.nzbin/2., 1., -1.)
			#	raise RuntimeError("sorry... increase args.nxbin or args.nzbin")
			# convolution and shift massf file
			convolve_data = ndimage.convolve(step_1d,massfrac_1d,mode='wrap')
			output = np.argmax(convolve_data)
			optimal_shift = int(len(step_1d)/2)-output
			#print(" optimal shift = {} at {} frame".format(optimal_shift,iframe))
			massfrac_1d_shift = np.roll(massfrac_1d, optimal_shift)
			totalmass_1d_shift = np.roll(totalmass_1d, optimal_shift)
			intrinsic_mf_profile[iframe] = intrinsic_mf_profile[iframe] + massfrac_1d_shift
			intrinsic_tm_profile[iframe] = intrinsic_tm_profile[iframe] + totalmass_1d_shift
	if n_vacant_times != 0:
		n_vacant_avg = n_vacant_avg + float(n_vacant_cells)/float(n_vacant_times)
		n_vacant_avg_times = n_vacant_avg_times + 1
	imod = hjung.time.process_print(iframe, n_frames, imod)
# print out stats
if n_vacant_times != 0:
	print("Once vacant cells occur,")
	print(" avg # vacant cells/frame = {}".format(float(n_vacant_avg)/float(n_vacant_avg_times)))
	#print(" avg # multi domain/frame = {}".format(n_multi_avg/float(n_multi_frame)))
	print("  (exclude frames with no multi domain case)") 

# delete frames which has continously vacant cells
#print(skip_frame_list)
if len(skip_frame_list) > 0:
	print("remove {} and {} frames which has continuously vacant cells or multi domains...".format(n_vacant_avg_times,len(skip_frame_list)-n_vacant_avg_times))
	skip_frame_list = np.array(skip_frame_list)
	intrinsic_mf_profile = np.delete(intrinsic_mf_profile,skip_frame_list,axis=0)
	intrinsic_tm_profile = np.delete(intrinsic_tm_profile,skip_frame_list,axis=0)

# average
intrinsic_mf_profile = intrinsic_mf_profile/float(args.nxbin*args.nxbin)
intrinsic_tm_profile = intrinsic_tm_profile/float(args.nxbin*args.nxbin)


## save number histogram trajectory
np.savetxt(args.omassf, intrinsic_mf_profile, 
	header='{} intrinsic mass1 fraction by molecules'.format(intrinsic_mf_profile.shape), fmt='%f', comments='# ')
np.save(args.omassf, intrinsic_mf_profile)
np.savetxt(args.otmass, intrinsic_tm_profile, 
	header='{}, total mass by molecules'.format(intrinsic_tm_profile.shape), fmt='%f', comments='# ')
np.save(args.otmass, intrinsic_tm_profile)

## timer
hjung.time.end_print(start_proc, start_prof)