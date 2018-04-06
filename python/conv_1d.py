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
# ver 2.0 - change program to give array for alignment and removing multilayers
#            and the alignemt and removing function are moved to align_remove_1d.py on 3/21/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Convolution for 1d mass fraction which provides optimal movement \
					and time indices to remove multilayer situation')
## args
parser.add_argument('-i', '--input', default='traj.massf', nargs='?', 
	help='raw mass fraction profile (npy file format) and exclude .npy in argument')
parser.add_argument('-crit', '--crit', default=0.0, nargs='?', type=float,
	help='criteria of autocorrelation to make a pulse function or get domain size [-1:1]')
parser.add_argument('-rm', '--remove', default='YES', nargs='?',
	help='Remove multi-layers trajectory? (YES/any)')
parser.add_argument('-half', '--half', default='YES', nargs='?',
	help='calculate domain size and save profiles of only last half or all? (YES/any)')
parser.add_argument('-file', '--file', default='YES', nargs='?',
	help='Do you need output files? Only if need domain size stats. please any except YES (YES/any)')
parser.add_argument('-o', '--output', default='.align', nargs='?', 
	help='output surfix for aligned profiles')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 2.0')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np

# default for args
args.omassf = args.input + args.output # save aligned mass fraction profiles
args.oacf   = args.input + '.acf'      # save autocorrelation function in 1D 
args.odsize = args.input + '.dsize'    # save domain sizes
args.oremov = args.input + '.remove'   # save iframes to remove
args.oopt   = args.input + '.opt'      # save optimal shift for alignment
args.input  = args.input + '.npy'
if abs(args.crit) > 1.0:
	raise ValueError(" wrong argument for args.crit {} but should be within [-1:1]".format(args.crit))

## timer
start_proc, start_prof = hjung.time.init()

## load data files
massfrac_1d_t = np.load(args.input)
nbin = len(massfrac_1d_t[0])
print(" #bin = {}".format(nbin))

## calculate autocorrelation function
acf_1d_t = hjung.analyze.autocorr_1d_t(massfrac_1d_t, 'wrap') 
slab_shift = int(len(acf_1d_t[0])/2.0)

## get a pulse function (sum of two step functions) from acf
criteria_acf = args.crit
step_1d_t = np.where(acf_1d_t - criteria_acf < 0.0, -1., 1.) # when we define domain size as zero points in acf and make pulse 

## find iframes with multilayers from pulse function
step_diff = np.diff(step_1d_t) # difference of neighbor element
n_frames = len(step_diff)
multilayer_iframes = [] 
for i_frame in range(n_frames):
	if 'YES' in args.half:
		if i_frame < int(n_frames/2):
			multilayer_iframes.append(i_frame)
			continue
	step_diff_iframe = step_diff[i_frame]
	step_n_up = int((step_diff_iframe > 0.0).sum())
	step_n_down = int((step_diff_iframe < 0.0).sum())
	# the sum of number of positives or negatives should have 1 if phase separation noramlly occurs
	if (step_n_up != 1) or (step_n_down != 1): 
		print(" multi-domain problem ({} {}) at {}".format(step_n_up,step_n_down,i_frame))
		multilayer_iframes.append(i_frame)
multilayer_iframes = np.array(multilayer_iframes,dtype=np.int)
print(" #frames to remove = {}".format(len(multilayer_iframes)))

## remove data using multilayer_iframes and args.half
massfrac_1d_t = np.delete(massfrac_1d_t,multilayer_iframes,axis=0)
acf_1d_t      = np.delete(acf_1d_t     ,multilayer_iframes,axis=0)
step_1d_t     = np.delete(step_1d_t    ,multilayer_iframes,axis=0)
step_diff     = np.delete(step_diff    ,multilayer_iframes,axis=0)	

# determine interface and domain size
step_up = np.argmax(step_diff,axis=1)
step_down = np.argmin(step_diff,axis=1)
domain_size = step_down - step_up
print("domain size (acf {}%) (avg,std) = {} +- {}".format(int(criteria_acf*100), np.average(domain_size), np.std(domain_size)))

## convolution and shift massf profile
if 'YES' in args.file:
	# convolution
	optimal_shift = hjung.analyze.convolve_1d_t(step_1d_t, massfrac_1d_t, 'wrap', 'max') 
	print(" optimal shift std = {}".format(np.std(optimal_shift)))
	n_frames = len(step_1d_t)
	align_mf = np.full_like(massfrac_1d_t,0.)
	# shifting
	for iframe in range(n_frames):
		align_mf[iframe] = np.roll(massfrac_1d_t[iframe], optimal_shift[iframe]) #optimal_shift[0]

	## write
	np.savetxt(args.oacf, acf_1d_t, 
		header='spatial autocorr(slab_lag,i_frame) (%d,%d) for delta_number, Plot u ($1-%d):2:3'	 
		%(len(acf_1d_t),len(acf_1d_t[0]),slab_shift), fmt='%f', comments='# ')
	np.savetxt(args.omassf, align_mf, 
		header='%d, %d, aligned massf fraction by ACF and molecules in nbins' \
		%(len(align_mf),nbin), fmt='%f', comments='# ')
	np.save(args.omassf, align_mf) 
	np.save(args.oremov, multilayer_iframes) 
	np.save(args.oopt, optimal_shift) 

## timer
hjung.time.end_print(start_proc, start_prof)