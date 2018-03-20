#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 2/17/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Get block average for input .npy file \n'
				'See details in Appendix D. Statistical Errors \n'
				' in the book Understanding Molecular Simulation \n'
				' by Daan Frenkel')
# args
parser.add_argument('-i', '--input', default='data.npy', nargs='?', 
	help='input file (only support 1d txt file, otherwise use .npy file)')
parser.add_argument('-o', '--output', default='.hist', nargs='?', 
	help='surfix of output file')
parser.add_argument('-b', '--begin', default=-1, nargs='?', type=int,
	help='beginning time iframe')
parser.add_argument('-e', '--end', default=-1, nargs='?', type=int,
	help='end time iframe')
parser.add_argument('-nbin', '--nbin', default=20, nargs='?', type=int,
	help='# bins')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import numpy as np
import math
import sys
sys.path.append('/home/htjung/Utility/python/')
import hjung
from hjung import *

## timer
start_proc, start_prof = hjung.time.init()

## load file
data = hjung.io.read_simple(args.input,args.begin,args.end)
n_frames = data.shape[0]
if len(data.shape) > 1:
	n_species = data.shape[1]
else:
	n_species = 1

## histogram per each species
hist_min = np.amin(data)
hist_max = np.amax(data)
data_hist = np.zeros((n_species,args.nbin))
i_species = 0
for i_data in np.transpose(data):
	data_hist[i_species], bin_edges = np.histogram(i_data,bins=args.nbin,range=(hist_min,hist_max),density=True)
	i_species += 1
data_hist = np.transpose(data_hist) 

# save
args.output = args.input.replace('.npy','') + args.output
data_set = np.column_stack((bin_edges[:-1],data_hist))
np.savetxt(args.output, data_set, 
	header='bin_edges, histogram per molecules', fmt='%e', comments='# ')
np.save(args.output, data_set)

# average
#print(np.transpose(data_hist)[0])
if n_species > 1:
	data_hist_avg = np.mean(data_hist,axis=1)
	data_hist_std = np.std(data_hist,axis=1)
	data_set = np.column_stack((bin_edges[:-1],data_hist_avg,data_hist_std))
	np.savetxt(args.output+str('.avg'), data_set,
		header='bin_edges, averaged histogram, std', fmt='%e', comments='# ')
	np.save(args.output+str('.avg'), data_set)

# make image
import matplotlib
matplotlib.use('Agg') # avoid the message "failed to get the current screen resources"
import matplotlib.pyplot as plt
width = np.diff(bin_edges)
center =  (bin_edges[:-1] + bin_edges[1:]) / 2
fig, ax = plt.subplots(figsize=(6,4))
if n_species > 1:
	ax.bar(center, data_hist_avg, align='center', width=width, yerr=data_hist_std)
else:
	ax.bar(center, np.transpose(data_hist)[0], align='center', width=width)
ax.axvline(x=np.mean(data),color='r',label='mean')
ax.axvline(x=np.median(data),color='b',label='median')
ax.legend()
#ax.set_xticks(bin_edges)
fig.savefig(args.output+str('.png'))

## timer
hjung.time.end_print(start_proc, start_prof)