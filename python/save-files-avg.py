#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 4/1/2018


import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='average files from .npy or text files which are the same size like (value, *)s' )
## args
parser.add_argument('-i', '--input', nargs='?', 
	help='input 1D profile file generated by .txt, .npy, else')
parser.add_argument('-num', '--number', default=0, nargs='?', type=int,
	help='index of files like args.input.[0:N-1]')
parser.add_argument('-o', '--output', default='INPUT', nargs='?', 
	help='output file of averaged 1D profile (.avg)')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.2')
# read args
args = parser.parse_args()
# check args
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np
import copy

# default for args
if 'INPUT' in args.output:
	new_name = args.input.replace('.npy','')
	new_name = new_name.replace('.txt','')
	args.output = new_name + '.avg'
else:
	args.output = args.output + '.avg'

## timer
start_proc, start_prof = hjung.time.init()

## read input file
for ifile in range(args.number):
	filename = args.input + '.' + str(ifile)
	if '.npy' in filename:
		data = np.load(filename)
	else:
		data = np.loadtxt(filename)
	#print(data.shape)
	#print(len(data))
	if 'data_collect' not in locals():
		data_collect = np.zeros((args.number,len(data))) # we only use numbers in the first column
	data_collect[ifile] = copy.copy(data[:,0])
data_mean = np.mean(data_collect,axis=0)
data_std = np.std(data_collect,axis=0)

data_prof = np.column_stack((data_mean,data_std))

## save averaged profile
print("="*30)
np.savetxt(args.output, data_prof, fmt='%f')
np.save(args.output, data_prof)
print("Finished saving averages.")

## timer
hjung.time.end_print(start_proc, start_prof)
