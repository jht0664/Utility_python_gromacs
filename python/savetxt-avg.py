#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on  2/26/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='block average 1D Profile from np.savetxt file')
## args
parser.add_argument('-i', '--input', nargs='?', 
	help='input 1D profile file generated by np.savetxt')
parser.add_argument('-b', '--begin', default=0, nargs='?', 
	help='index of beginning frame [0:N-1]')
parser.add_argument('-e', '--end', default=-1, nargs='?', 
	help='index of end frame [0:N-1]. If negative, use end frame')
parser.add_argument('-tol', '--tol', default=0.0, nargs='?', 
	help='tolerance for block average (> 0 and <= 1). No block average (tol=0). # frames to average (tol>1)')
parser.add_argument('-o', '--output', default='INPUT.avg', nargs='?', 
	help='output file of block averaged 1D profile')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# default for args
if args.output is 'INPUT.avg':
	args.output = args.input+'.avg'
args.begin = int(args.begin)
args.end = int(args.end)
args.tol = float(args.tol)

## Check arguments for log
print("===============================")
print("input filename   = ", args.input)
print("index of beginning frame = ", args.begin)
if args.end != -1:
	print("index of end frame = ", args.end)
else:
	print("Set end frame is the end.")
if args.tol == 0.0:
	print("Set no bloack average")
elif args.tol <= 1.0: 
	print("tolerance for block average = %f" %args.tol)
elif args.tol > 1.0:
	print("set block length = %d" %(int(args.tol)))
print("output filename = ", args.output)

## timer
import time
start_clock = time.clock() # process time
start_wall = time.time() # wall time

## import modules
import hjung
from hjung import *
import numpy as np

## check argument
if args.tol < 0.0:
	raise ValueError("wrong input of tolerance, %f" %args.tol)
elif args.tol >= 1.0: 
	print("Warning: tolerance %f is assigned to block_size, %d" %(args.tol, int(args.tol)))
else:
	print("="*30)

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

## block average
print("="*30)
block_length = hjung.analyze.opt_block_length_1d_t(data_1d_t,args.tol) 
data_1d_t = hjung.analyze.block_average_1d(data_1d_t,block_length)
print("Done: Block average")

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

print("="*30)
print(time.clock() - start_clock, "seconds process time")
print(time.time() - start_wall, "seconds wall time")