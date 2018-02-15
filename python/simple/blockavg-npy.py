#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 9/08/2016
#         - (Ref.) Appendix D. Statistical Errors 
#       		in the book Understanding Molecular Simulation 
# 				by Daan Frenkel
#		  - (Ref.) Chapter2 in annual Reports in Computational Chemistry, 
#				Vol 5, 2009, doi: 10.1016/S1574-1400(09)00502-7

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Get block average for input .npy file \n'
				'See details in Appendix D. Statistical Errors \n'
				' in the book Understanding Molecular Simulation \n'
				' by Daan Frenkel')
# args
parser.add_argument('-i', '--input', default='block.npy', nargs='?', 
	help='input file')
parser.add_argument('-o', '--output', default='.bavg', nargs='?', 
	help='surfix of output file')
parser.add_argument('-b', '--begin', default=-1, nargs='?', type=int,
	help='beginning time iframe')
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
data = np.load(args.input)
n_frames = data.shape[0]
if args.begin == -1:
	data = data[int(n_frames/2):]
else:
	data = data[args.begin:]
n_frames = data.shape[0]
n_species = data.shape[1] # number of independent molecules to average
print("total {} frames and {} species".format(n_frames,n_species))

# make a list of block length we will do
list_block_length = []
block_length = 5
while True:
	num_block = int(n_frames/block_length)
	if num_block > 5:
		if n_frames%block_length == 0:
			list_block_length.append(block_length)
		block_length = block_length + 1
		continue
	else:
		break
list_block_length = np.array(list_block_length)

# computation of block length, # blocks, blocked average,
#  and blocked standard error (= square root of variance = standard deviation)
data_avg = np.mean(data,axis=0)
avg = np.mean(data_avg)
data_std = np.std(data,axis=0)
std = np.std(data_avg)

list_bse = []
for block_length in list_block_length:
	list_bavg = []
	for i in range(0,len(data),block_length):
		# mean of each block
		bavg_list = np.mean(data[i:i+block_length-1],axis=0)
		bavg = np.mean(bavg_list)
		list_bavg.append(bavg)
	list_bavg = np.array(list_bavg)
	# assume all blocks have the same weight
	bstd = np.std(list_bavg)
	bavg_b = np.mean(list_bavg)
	bse = bstd/math.sqrt(len(list_bavg))
	list_bse.append([block_length, bse, bavg_b, bstd])
list_bse = np.array(list_bse)


np.savetxt(args.output, list_bse, 
	header='# block length, blocked standard error, blocked average, blocked STD (avg = {} +- {}'.format(avg,std), fmt='%f', comments='# ')
np.save(args.output, list_bse)

## timer
hjung.time.end_print(start_proc, start_prof)