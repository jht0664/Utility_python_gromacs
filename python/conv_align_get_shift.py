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
	description='get convolution alignment shift array using a profile')
## args
parser.add_argument('-i', '--input', default='traj.mass', nargs='?', 
	help='raw mass profile (npy file format) and exclude .npy in argument')
parser.add_argument('-rm', '--remove', default='YES', nargs='?',
	help='Remove multi-layers trajectory? (YES/any)')
parser.add_argument('-half', '--half', default='YES', nargs='?',
	help='calculate domain size and save profiles of only last half or all? (YES/any)')
parser.add_argument('-o', '--output', default='.conv', nargs='?', 
	help='output surfix for shift array')
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

# default for args
args.oconv = args.input + args.output # save shift array (how bins it should shift for alignment)
args.oacf = args.input + '.acfs'       # save autocorrelation function
args.oremoveframe = args.input + '.removeframe'   # save iframes when multilayers happen 
args.input = args.input + '.npy'

## timer
start_proc, start_prof = hjung.time.init()

## load data files
mass_1d_t = np.load(args.input)
nbin = len(mass_1d_t[0])

## calculate autocorrelation function
acf_1d_t_wrap = hjung.analyze.autocorr_1d_t(mass_1d_t, 'wrap') 
slab_shift = int(len(acf_1d_t_wrap[0])/2.0)
np.savetxt(args.oacf, acf_1d_t_wrap, 
	header='spatial autocorr(slab_lag,i_frame) (%d,%d) for delta_number, Plot u ($1-%d):2:3'	 
	%(len(acf_1d_t_wrap),len(acf_1d_t_wrap[0]),slab_shift), fmt='%f', comments='# ')

## get shift array
step_1d_t_wrap = np.where(acf_1d_t_wrap < 0.0, -1., 1.) # when we define domain size as zero points in acf
align_shift = hjung.analyze.convolve_1d_t(step_1d_t_wrap, mass_1d_t, 'wrap', 'max') 

## multilayer check
def multilayer_in_step_fn(step_1d_t_wrap):
	step_diff = np.diff(step_1d_t_wrap) # difference of neighbor element
	# save iframes when multilayer occurs
	multilayer_iframes = [] 
	for i_frame in range(len(step_diff)):
		step_diff_iframe = step_diff[i_frame]
		step_n_up = (step_diff_iframe > 0.0).sum() 
		step_n_down = (step_diff_iframe < 0.0).sum()
		if (step_n_up > 1.0) or (step_n_down > 1.0):
			#print("Probably it has two or more layers (multi-domains) at {}. We remove them in profiles!".format(i_frame))
			multilayer_iframes.append(i_frame)
	multilayer_iframes = np.array(multilayer_iframes,dtype=np.int)
	# determine interface 
	step_up = np.argmax(step_diff,axis=1)
	step_down = np.argmin(step_diff,axis=1)
	return step_up, step_down, multilayer_iframes

def remove_data(step_up, step_down, multilayer_iframes, ask_half):
	# remove all of multilayers
	step_up = np.delete(step_up,multilayer_iframes,axis=0)	
	step_down = np.delete(step_down,multilayer_iframes,axis=0)	
	# remove first half trajectories
	if 'YES' in ask_half:
		remove_range = np.arange(len(step_up)/2)
		step_up = np.delete(step_up,remove_range,axis=0)	
		step_down = np.delete(step_down,remove_range,axis=0)	
	return step_up, step_down

def print_domain_size(step_up, step_down, text1):
	domain_size = step_down - step_up
	domain_size_avg = np.average(domain_size)
	domain_size_std = np.std(domain_size)
	print("domain size {} (avg,std) = {} +- {}".format(text1, domain_size_avg, domain_size_std))
	return domain_size

def main_domain_size_step_fn(acf_1d_t, criteria_massf, ask_half, text_print):  
	step_1d_t = np.where(acf_1d_t - criteria_massf < 0.0, -1., 1.) # when we define domain size as zero points in acf
	step_up, step_down, multilayer_iframes = multilayer_in_step_fn(step_1d_t)
	step_up, step_down = remove_data(step_up, step_down, multilayer_iframes, ask_half)
	domain_size = print_domain_size(step_up, step_down, text_print)
	return domain_size, multilayer_iframes

domain_size, multilayer_iframes = main_domain_size_step_fn(acf_1d_t_wrap, 0.5, args.half, "(50%)")
domain_size, multilayer_iframes = main_domain_size_step_fn(acf_1d_t_wrap, 0.,  args.half, "(0%)")
print(" removed {} frames ({:.3}%) due to multilayers".format(len(multilayer_iframes),len(multilayer_iframes)*100/len(acf_1d_t_wrap)))
if len(multilayer_iframes) > int(len(acf_1d_t_wrap)/50):
	raise RuntimeError("# frames to have multilayers > half frames of trajectory. Check your trajectory.")

## remove all of multilayers
## remove first half trajectories
removeframe = []
removeframe = np.array(removeframe,dtype=int)
if 'YES' in args.remove:
	removeframe = np.append(removeframe, multilayer_iframes)
if 'YES' in args.half:
	removeframe = np.append(removeframe, np.arange(len(acf_1d_t_wrap)/2))

## write
np.savetxt(args.oconv, align_shift, 
	header='shift array in unit of bins', fmt='%f', comments='# ')
np.save(args.oconv, align_shift) 
np.savetxt(args.oremoveframe, removeframe, 
	header='iframe list to remove due to multilayers and constraints', fmt='%f', comments='# ')
np.save(args.oremoveframe, removeframe) 

## timer
hjung.time.end_print(start_proc, start_prof)