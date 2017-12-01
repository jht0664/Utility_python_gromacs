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

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Convolution alignment for 1d mass fraction and total mass profiles')
## args
parser.add_argument('-imf', '--in_massf', default='traj.massf', nargs='?', 
	help='raw mass fraction profile (npy file format) and exclude .npy in argument')
parser.add_argument('-itm', '--in_tmass', default='traj.tmass', nargs='?', 
	help='raw totmal mass profile (npy file format) and exclude .npy in argument')
parser.add_argument('-rm', '--remove', default='YES', nargs='?',
	help='Remove multi-layers trajectory? (YES/any)')
parser.add_argument('-half', '--half', default='YES', nargs='?',
	help='calculate domain size and save profiles of only last half or all? (YES/any)')
parser.add_argument('-o', '--output', default='.align', nargs='?', 
	help='output surfix for aligned profiles')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.2')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np

# default for args
args.omassf = args.in_massf + '.align' # save aligned mass fraction profiles
args.otmass = args.in_tmass + '.align' # save aligned total mass profiles
args.oacf = args.in_massf + '.acf'     # save autocorrelation function in 1D 
args.odsize = args.in_massf + '.dsize' # save domain sizes
args.omulti = args.in_massf + '.multi' # save iframes when multilayers happen 
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

## convert to autocorrelation function to step-function
step_1d_t_wrap = np.where(acf_1d_t_wrap < 0.0, -1., 1.)
align_massfrac_1d_t, align_totalmass_1d_t = hjung.analyze.align_acf_w_data2(massfrac_1d_t, totalmass_1d_t, step_1d_t_wrap, 'wrap') 

## find multi layers and remove all
step_diff = np.diff(step_1d_t_wrap)
multilayer_iframes = []
for i_frame in range(len(step_diff)):
	step_diff_iframe = step_diff[i_frame]
	step_n_up = (step_diff_iframe > 0.0).sum()
	step_n_down = (step_diff_iframe < 0.0).sum()
	if (step_n_up > 1.0) or (step_n_down > 1.0):
		print("Probably it has two or more layers (multi-domains) at {}. We remove them in profiles!".format(i_frame))
		multilayer_iframes.append(i_frame)
multilayer_iframes = np.array(multilayer_iframes,dtype=np.int)
step_up = np.argmax(step_diff,axis=1)
step_down = np.argmin(step_diff,axis=1)
# remove all of multilayers
if 'YES' in args.remove:
	step_up = np.delete(step_up,multilayer_iframes,axis=0)	
	step_down = np.delete(step_down,multilayer_iframes,axis=0)	
	align_massfrac_1d_t = np.delete(align_massfrac_1d_t,multilayer_iframes,axis=0)	
	align_totalmass_1d_t = np.delete(align_totalmass_1d_t,multilayer_iframes,axis=0)	

## domain size information
# remove first half trajectories
if 'YES' in args.half:
	remove_range = np.arange(len(step_up)/2)
	step_up = np.delete(step_up,remove_range,axis=0)	
	step_down = np.delete(step_down,remove_range,axis=0)	
	align_massfrac_1d_t = np.delete(align_massfrac_1d_t,remove_range,axis=0)	
	align_totalmass_1d_t = np.delete(align_totalmass_1d_t,remove_range,axis=0)	
domain_size = step_down - step_up
domain_size_avg = np.average(domain_size)
domain_size_std = np.std(domain_size)
print("domain size (avg,std) = {0:.5f} +- {1:.5f}".format(domain_size_avg,domain_size_std))

## write
np.savetxt(args.omassf, align_massfrac_1d_t, 
	header='%d, %d, aligned mass1 fraction by ACF and molecules in nbins' \
	%(len(align_massfrac_1d_t),nbin), fmt='%f', comments='# ')
#np.save(args.omassf, align_massfrac_1d_t) 
np.savetxt(args.otmass, align_totalmass_1d_t, 
	header='%d, %d, aligned total mass by ACF and molecules in nbins' \
	%(len(align_totalmass_1d_t),nbin), fmt='%f', comments='# ')
#np.save(args.otmass, align_totalmass_1d_t) 
np.savetxt(args.odsize, domain_size,
	header='domain size in #bins from step-wise ACF, (avg,std) = %d +- %d' \
	%(domain_size_avg,domain_size_std), fmt='%f', comments='# ')
#np.save(args.odsize, domain_size) 
np.savetxt(args.omulti, multilayer_iframes,
	header='iframes when multilayers occurs', fmt='%d \n', comments='# ')

## timer
hjung.time.end_print(start_proc, start_prof)