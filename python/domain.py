#!/usr/bin/env python3
# ver 0.1 - copy from conv_align.py (v1.3) and remove alignment function on 12/10/2017
# ver 0.2 - change the purpose of this program and modify on 01/21/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='get domain info. using convolution of a given (mass) profile')
## args
parser.add_argument('-i', '--input', default='traj.mass', nargs='?', 
	help='raw mass profile (npy file format) and exclude .npy in argument')
parser.add_argument('-b', '--begin', default=-1, nargs='?', type=int,
	help='beginning index of frames (-1 means index of middle number of frames) ')
parser.add_argument('-step', '--step', default=1, nargs='?', type=int,
	help='step index of frames')
parser.add_argument('-e', '--end', default=-1, nargs='?', type=int,
	help='end index of frames (-1 means the index of last frames) ')
parser.add_argument('-max_frame', '--max_frame', default=50, nargs='?', type=int,
	help='collect up to the maximum frames for each rdfs ')
parser.add_argument('-o', '--output', default='.domain', nargs='?', 
	help='output surfix for dataset [iframes, check_multilayer(Y:1,N:0), align_shift, step_up, step_down, domain_size]')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.2')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np
import random as rd
import itertools

# default for args
args.oset = args.input + args.output  # save [iframe/domain size/interface position]
args.odic = args.input + args.output + '.dic' # save dictionary
args.input = args.input + '.npy'
if args.begin > args.end:
	raise ValueError("args.begin is wrong.")

## timer
start_proc, start_prof = hjung.time.init()

## load data files
mass_1d_t = np.load(args.input)
nbin = len(mass_1d_t[0])

## calculate autocorrelation function
acf_1d_t_wrap = hjung.analyze.autocorr_1d_t(mass_1d_t, 'wrap') 
slab_shift = int(len(acf_1d_t_wrap[0])/2.0)
#np.savetxt(args.oacf, acf_1d_t_wrap, 
#	header='spatial autocorr(slab_lag,i_frame) (%d,%d) for delta_number, Plot u ($1-%d):2:3'	 
#	%(len(acf_1d_t_wrap),len(acf_1d_t_wrap[0]),slab_shift), fmt='%f', comments='# ')

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

def main_domain_size_step_fn(acf_1d_t, criteria_massf):  
	step_1d_t = np.where(acf_1d_t - criteria_massf < 0.0, -1., 1.) # when we define domain size as zero points in acf
	step_up, step_down, multilayer_iframes = multilayer_in_step_fn(step_1d_t) # get iframes when multilayers happen, and step_up and step_down positions
	return step_up, step_down, multilayer_iframes

step_up, step_down, multilayer_iframes = main_domain_size_step_fn(acf_1d_t_wrap, 0.5) # set criteria acf 50% 
print(" multilayer {} frames ({:.3}%)".format(len(multilayer_iframes),len(multilayer_iframes)*100/len(acf_1d_t_wrap)))

check_multilayer = np.zeros(np.size(step_up),dtype=int)
for iframe in multilayer_iframes:
	check_multilayer[iframe] = 1 # if multilayer occurs at iframe, check_multilayer save the value 1. Otherwise, 0.
domain_size = step_down - step_up

n_frames = np.size(step_up)
iframes = np.arange(np.size(step_up),dtype=int)

## select frames to collect data
if args.end > n_frames:
	raise ValueError("args.end is greater than total number of frames")
if args.end == -1:
	args.end = n_frames
if args.begin == -1:
	args.begin = int(args.end/2.0)
elif args.begin < -1:
	args.begin = int(args.end + args.begin)
if args.step < 1:
	raise ValueError("args.step is wrong")
print("set possible frame range = [{}, {}, {}]".format(args.begin,args.end,args.step))
possible_n_frames = np.arange(start=args.begin,stop=args.end,step=args.step)

## reassign arrays to reduce memory usage
iframes_list = np.setdiff1d(possible_n_frames,np.nonzero(check_multilayer)) # get iframe list without multilayer frames
print("reduced {} frames due to multilayers".format(len(possible_n_frames)-len(iframes_list)))
n_frames = len(iframes_list)
print("Finally, save data for total {} frames".format(n_frames))

## check domain size distrbituion
domain_size_stat = domain_size[iframes_list]
domain_size_median = int(np.median(domain_size_stat))
print("Domain size = {} (median)".format(domain_size_median))
domain_list, domain_counts = np.unique(domain_size_stat, return_counts=True)
print("Domain size counts = {}".format(dict(zip(domain_list, domain_counts))))
domain_dens = domain_counts/len(domain_size_stat)
print("Domain size distribution = {}".format(dict(zip(domain_list, domain_dens))))
# draw histogram of domain_size
import matplotlib.pyplot as plt
plt.switch_backend('agg')
# the histogram of the data
plt.figure()
n, bins, patches = plt.hist(domain_size_stat, int(len(domain_list)+1), color='green')
print(bins)
# add a 'best fit' line
plt.xlabel('Domain sizes')
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ Domain Sizes}\ $')
plt.grid(True)
plt.savefig(args.odic+'.png')
## way2: determine number of partial rdfs and check symmetry of domain sizes
n_domain_sizes = len(domain_list)
i_rdf = np.zeros(n_domain_sizes,dtype=int)
if n_domain_sizes > 4:
	avg_domain = np.mean(domain_size_stat)
	for i_domain in range(n_domain_sizes):
		if domain_dens[i_domain] < 0.25:
			if domain_list[i_domain] < avg_domain:
				i_rdf[i_domain] = 1
			elif domain_list[i_domain] > avg_domain:
				i_rdf[i_domain] = 2
			else:
				raise ValueError(" something wrong to decide i_rdf")
## way1: determine number of partial rdfs and check symmetry of domain sizes
#if len(domain_list)%2 == 1: # odd number of domain sizes
#	if domain_list[int((len(domain_list)-1)/2)] != domain_size_median:
#		print(" exceptional case: median of (odd) domain size list is not at center")
#	# make partial rdfs with every 2 bins far from median 
#	n_rdfs = int((len(domain_list)+1)/2)
#else: # even number of domain sizes
#	if len(domain_list) == 2:
#		raise ValueError(" too small pool of domain sizes.. need to run longer simulations")
#	if (domain_list[int(len(domain_list)/2)-1] != domain_size_median) and \
#		 (domain_list[int(len(domain_list)/2)] != domain_size_median):
#		print(" exceptional case: median of (even) domain size list is not near center")
#	n_rdfs = int(len(domain_list)/2)
#	domain_size_median = (domain_list[int(len(domain_list)/2)-1] + domain_list[int(len(domain_list)/2)])/2.0
#	print("set a new median domain size  = {}".format(domain_size_median))
#
#i_rdf = np.floor(np.absolute(domain_list - domain_size_median)/2.0)
dict_domain = dict(zip(domain_list, i_rdf))
print("Here is dict: {}".format(dict_domain))

## limit data within args.max_frame
select_frames = []
for i_domain_size in domain_list:
	#print(domain_size_stat,i_domain_size)
	specific_domain_list = np.where(domain_size_stat == i_domain_size)
	#print(specific_domain_list[0])
	#print(len(specific_domain_list[0]))
	over_frames = len(specific_domain_list[0]) - args.max_frame
	#print(over_frames)
	if over_frames > 0:
		select_tmp = rd.sample(list(specific_domain_list[0]),over_frames)
		select_frames.append(select_tmp)
		#print(len(select_tmp),len(select_frames))
select_frames = list(itertools.chain.from_iterable(select_frames)) # flatten the list of list
#print(len(select_frames))
iframes_list = np.delete(iframes_list, select_frames)
dataset = np.column_stack((iframes_list,check_multilayer[iframes_list],
	align_shift[iframes_list],step_up[iframes_list],step_down[iframes_list],domain_size[iframes_list]))

## write
#np.savetxt(args.odic, dict_domain, fmt='%i')
with open(args.odic,'w') as file_dic:
	file_dic.write(str(dict_domain))
np.savetxt(args.oset, dataset, 
	header=' [iframes, check_multilayer(Y:1,N:0), align_shift, step_up, step_down, domain_size] in unit of bins', fmt='%d', comments='# ')
np.save(args.oset, dataset) 

## timer
hjung.time.end_print(start_proc, start_prof)