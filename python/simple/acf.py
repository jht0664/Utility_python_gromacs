#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 2/17/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='autocorrelation functions for each columns (molecules)')
# args
parser.add_argument('-i', '--input', default='acf.npy', nargs='?', 
	help='input .npy file')
parser.add_argument('-o', '--output', default='.acf', nargs='?', 
	help='surfix of output file')
parser.add_argument('-b', '--begin', default=-1, nargs='?', type=int,
	help='beginning time iframe')
parser.add_argument('-e', '--end', default=-1, nargs='?', type=int,
	help='end time iframe')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

args.output = args.input + args.output

## import modules
import numpy as np
import math
import sys
sys.path.append('/home/htjung/Utility/python/')
import hjung
from hjung import *

## timer
start_proc, start_prof = hjung.time.init()

## load rdf files
data = hjung.io.read_simple(args.input,args.begin,args.end)
n_frames = data.shape[0]
if len(data.shape) > 1:
	n_species = data.shape[1]
else:
	n_species = 1


# autocorrelation function
acf_data = hjung.analyze.autocorr_1d_t(np.transpose(data),'constant') # assume no correlation between start- and end-points
acf_out = np.transpose(acf_data)
acf_out = acf_out[int(n_frames/2):] # remove half due to symmetry

#acf_data = hjung.analyze.autocorr_signal_1d_t(np.transpose(data),'full')
#acf_out = np.transpose(acf_data)
#acf_out2 = acf_out2[int(n_frames/2):] # remove half due to symmetry
#print(acf_out2.shape)
#print(acf_out2)

# save data and avg
args.output = args.input.replace('.npy','') + args.output
np.savetxt(args.output, acf_out, 
	header='acfs for {} species with {} frames'.format(n_species,n_frames), fmt='%f', comments='# ')
np.save(args.output, acf_out)
acf_out_avg = np.mean(acf_out,axis=1)
acf_out_std = np.std(acf_out,axis=1)
avg_data = np.column_stack((acf_out_avg,acf_out_std))
np.savetxt(args.output+str('.avg'), avg_data, 
	header='averaged acf (mean, std)', fmt='%f', comments='# ')
np.save(args.output+str('.avg'), avg_data)


## timer
hjung.time.end_print(start_proc, start_prof)
