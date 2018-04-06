#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 03/21/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Convolution alignment for 1d mass fraction and total mass profiles')
## args
parser.add_argument('-i', '--input', default='traj.tmass', nargs='?', 
	help='raw profile (npy file format)')
parser.add_argument('-irem', '--irem', default='traj.massf.remove', nargs='?', 
	help='array of iframes to remove (npy file format)')
parser.add_argument('-iopt', '--iopt', default='traj.massf.opt', nargs='?', 
	help='array of optimal shift from conv_1d.py')
parser.add_argument('-o', '--output', default='.align', nargs='?', 
	help='output surfix for aligned profiles')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np

# default for args
args.output = args.input + args.output # save aligned total mass profiles
args.input  = args.input + '.npy'
args.irem   = args.irem + '.npy'
args.iopt   = args.iopt + '.npy'

## timer
start_proc, start_prof = hjung.time.init()

## load data files
data_1d_t = np.load(args.input)
multilayer_iframes = np.load(args.irem)
optimal_shift = np.load(args.iopt)

## remove time frames
data_1d_t = np.delete(data_1d_t,multilayer_iframes,axis=0)
if len(optimal_shift) != len(data_1d_t):
	raise ValueError(" different array sizes {} {}".format(len(optimal_shift),len(data_1d_t)))

## move by optimal shift
n_frames = len(data_1d_t)
align_data = np.full_like(data_1d_t,0.)
for iframe in range(n_frames):
	align_data[iframe] = np.roll(data_1d_t[iframe], optimal_shift[iframe]) #optimal_shift[0]

## write
np.savetxt(args.output, align_data, 
	header='%d, %d, aligned profile by optimal shift' \
	%(len(align_data),len(align_data[0])), fmt='%f', comments='# ')
np.save(args.output, align_data) 

## timer
hjung.time.end_print(start_proc, start_prof)