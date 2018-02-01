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

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Convolution alignment for 1d mass fraction and total mass profiles')
## args
parser.add_argument('-imf', '--in_massf', default='traj.massf', nargs='?', 
	help='raw mass fraction profile (npy file format) and exclude .npy in argument')
parser.add_argument('-itm', '--in_tmass', default='traj.tmass', nargs='?', 
	help='raw totmal mass or mass profile (npy file format) and exclude .npy in argument')
parser.add_argument('-rm', '--remove', default='YES', nargs='?',
	help='Remove multi-layers trajectory? (YES/any)')
parser.add_argument('-half', '--half', default='YES', nargs='?',
	help='calculate domain size and save profiles of only last half or all? (YES/any)')
parser.add_argument('-o', '--output', default='.align', nargs='?', 
	help='output surfix for aligned profiles')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.3')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np
import matplotlib
matplotlib.use('Agg') # avoid to show figures when running bash shell script
import matplotlib.pyplot as plt
from scipy import stats
import copy

# default for args
args.omassf = args.in_massf + args.output # save aligned mass fraction profiles
args.otmass = args.in_tmass + args.output # save aligned total mass profiles
args.oacf = args.in_massf + '.acf'        # save autocorrelation function in 1D 
args.odsize = args.in_massf + '.dsize'    # save domain sizes
#args.omulti = args.in_massf + '.multi'    # save iframes when multilayers happen 
args.in_massf = args.in_massf + '.npy'
args.in_tmass = args.in_tmass + '.npy'

## timer
start_proc, start_prof = hjung.time.init()

## load data files
massfrac_1d_t = np.load(args.in_massf)
totalmass_1d_t = np.load(args.in_tmass)
if massfrac_1d_t.size != totalmass_1d_t.size:
	raise ValueError("the size of two data files are different.")
nbin = len(massfrac_1d_t[0])

## calculate autocorrelation function
acf_1d_t_wrap = hjung.analyze.autocorr_1d_t(massfrac_1d_t, 'wrap') 
slab_shift = int(len(acf_1d_t_wrap[0])/2.0)
np.savetxt(args.oacf, acf_1d_t_wrap, 
	header='spatial autocorr(slab_lag,i_frame) (%d,%d) for delta_number, Plot u ($1-%d):2:3'	 
	%(len(acf_1d_t_wrap),len(acf_1d_t_wrap[0]),slab_shift), fmt='%f', comments='# ')

## convert to autocorrelation function to step-function and align density profiles
step_1d_t_wrap = np.where(acf_1d_t_wrap < 0.0, -1., 1.) # when we define domain size as zero points in acf
align_massfrac_1d_t, align_totalmass_1d_t = hjung.analyze.align_acf_w_data2(massfrac_1d_t, totalmass_1d_t, step_1d_t_wrap, 'wrap') 

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
			print("Probably it has two or more layers (multi-domains) at {}. We remove them in profiles!".format(i_frame))
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

domain_size, multilayer_iframes = main_domain_size_step_fn(acf_1d_t_wrap, 0.2, args.half, "(20%)")
domain_size, multilayer_iframes = main_domain_size_step_fn(acf_1d_t_wrap, 0.8, args.half, "(80%)")
domain_size, multilayer_iframes = main_domain_size_step_fn(acf_1d_t_wrap, 0.5, args.half, "(50%)")

## relationship between domain size (50%), massf, and tmass 
# however, when I checked the result, no correlation because of Monte Carlo simulation. 
# (only average is meaningful and one movement does not affect rest of part of system, i.e. environment)
center = int(len(align_massfrac_1d_t[0])/2 - 1)
massf_center = copy.copy(align_massfrac_1d_t[:,center])
tmass_center = copy.copy(align_totalmass_1d_t[:,center])
massf_center, tmass_center = remove_data(massf_center,tmass_center,multilayer_iframes,args.half)
## domain size - massf (not necessary jobs)
#plt.figure()
#dm_s, dm_i, dm_r, dm_p, dm_err = stats.linregress(domain_size, massf_center)
#print("r-squared (domain-massf) = {}".format(dm_r**2.))
#plt.plot(domain_size, massf_center, 'o', label='data')
#plt.plot(domain_size, dm_i + dm_s*domain_size, 'r', label='fit')
#plt.legend()
#plt.savefig(args.odsize+'.dm.png')
## domain size - tmass
#plt.figure()
#dm_s, dm_i, dm_r, dm_p, dm_err = stats.linregress(domain_size, tmass_center)
#print("r-squared (domain-tmass) = {}".format(dm_r**2.))
#plt.plot(domain_size, tmass_center, 'o', label='data')
#plt.plot(domain_size, dm_i + dm_s*domain_size, 'r', label='fit')
#plt.legend()
#plt.savefig(args.odsize+'.dt.png')
## tmass - massf
#plt.figure()
#dm_s, dm_i, dm_r, dm_p, dm_err = stats.linregress(tmass_center, massf_center)
#print("r-squared (tmass-massf) = {}".format(dm_r**2.))
#plt.plot(tmass_center, massf_center, 'o', label='data')
#plt.plot(tmass_center, dm_i + dm_s*tmass_center, 'r', label='fit')
#plt.legend()
#plt.savefig(args.odsize+'.tm.png')
## save array stacks for output
domainsize_massf_tmass = np.column_stack((domain_size, massf_center, tmass_center))

## remove all of multilayers
domain_size, multilayer_iframes = main_domain_size_step_fn(acf_1d_t_wrap, 0.,  args.half, "(0%)")
if 'YES' in args.remove:
	align_massfrac_1d_t = np.delete(align_massfrac_1d_t,multilayer_iframes,axis=0)	
	align_totalmass_1d_t = np.delete(align_totalmass_1d_t,multilayer_iframes,axis=0)	
## remove first half trajectories
if 'YES' in args.half:
	remove_range = np.arange(len(align_massfrac_1d_t)/2)
	align_massfrac_1d_t = np.delete(align_massfrac_1d_t,remove_range,axis=0)	
	align_totalmass_1d_t = np.delete(align_totalmass_1d_t,remove_range,axis=0)	

## write
np.savetxt(args.omassf, align_massfrac_1d_t, 
	header='%d, %d, aligned massf fraction by ACF and molecules in nbins' \
	%(len(align_massfrac_1d_t),nbin), fmt='%f', comments='# ')
#np.save(args.omassf, align_massfrac_1d_t) 
np.savetxt(args.otmass, align_totalmass_1d_t, 
	header='%d, %d, aligned (total or selected) mass by ACF and molecules in nbins' \
	%(len(align_totalmass_1d_t),nbin), fmt='%f', comments='# ')
#np.save(args.otmass, align_totalmass_1d_t) 
np.savetxt(args.odsize, domainsize_massf_tmass,
	header='domain size, mass fraction and total mass in center of domain', fmt='%f', comments='# ')
#np.save(args.odsize, domainsize_massf_tmass) 
#np.savetxt(args.omulti, multilayer_iframes,
#	header='iframes when multilayers occurs', fmt='%d \n', comments='# ')

## timer
hjung.time.end_print(start_proc, start_prof)