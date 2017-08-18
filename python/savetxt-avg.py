#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on  2/26/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='block average 1D Profile from np.savetxt file')
## args
parser.add_argument('-i', '--input', nargs='?', 
	help='input 1D profile file generated by np.savetxt')
parser.add_argument('-b', '--begin', default=0, nargs='?', type=int,
	help='index of beginning frame [0:N-1]')
parser.add_argument('-e', '--end', default=-1, nargs='?', type=int,
	help='index of end frame [0:N-1]. If negative, use end frame')
parser.add_argument('-tol', '--tol', default=0.0, nargs='?', type=float,
	help='tolerance for block average (> 0 and <= 1). No block average (tol=0). # frames to average (tol>1)')
parser.add_argument('-o', '--output', default='INPUT.avg', nargs='?', 
	help='output file of block averaged 1D profile')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()

## import modules
import hjung
from hjung import *
import numpy as np

# default for args
if args.output is 'INPUT.avg':
	args.output = args.input+'.avg'
args.begin = int(args.begin)
args.end = int(args.end)

## Check arguments for log
print("===============================")
print("input filename   = ", args.input)
print("index of beginning frame = ", args.begin)
if args.end != -1:
	print("index of end frame = ", args.end)
else:
	print("Set end frame is the end.")
hjung.blockavg.print_init(args.tol)
print("output filename = ", args.output)

## timer
start_proc, start_prof = hjung.time.init()

## check argument
args.tol = hjung.blockavg.check(args.tol)

## read input file
print("="*30)
data_1d_t = np.loadtxt(args.input, comments='#')
print("Total trajecotry has %d frames." %(len(data_1d_t)))
if args.end >= len(data_1d_t):
	raise ValueError("end frame is beyond trajectory")
if args.end < 0:
	data_1d_t = data_1d_t[args.begin:]
else:
	data_1d_t = data_1d_t[args.begin:args.end]
print("Done: reading input file")

## block average to get stable volume (or number density)
print("="*30)
data_1d_t, block_length = hjung.blockavg.main_1d(data_1d_t, None, args.tol) 

## make average and std
# use numpy.mean(array,axis=0) to avoid transpose cost
data_1d_avg = np.mean(data_1d_t, axis=0)
data_1d_std = np.std(data_1d_t, axis=0)
data_1d = np.column_stack((data_1d_avg,data_1d_std))

## save averaged profile
print("="*30)
np.savetxt(args.output, data_1d, 
	header='begin frame = %d, end frame = %d, generated by Hyuntae python code' %(args.begin,args.end), fmt='%f', comments='# ')
print("Finished saving average number density.")

## timer
hjung.time.end_print(start_proc, start_prof)