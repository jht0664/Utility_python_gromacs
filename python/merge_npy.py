#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 04/06/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='merge many .npy files (2d array) into one numpy file')
## args
parser.add_argument('-o', '--output', default='merge_file', nargs='?', 
	help='output prefix')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))


import glob
import os
import copy
import numpy as np

start = 0
npyfiles = glob.glob("*.npy")
npyfiles.sort()
print(npyfiles)
n_files = len(npyfiles)
print(n_files)
if n_files < 1:
	raise ValueError("no .npy files")
for npyfile in npyfiles:
	# get path
	filepath = os.path.join('./', npyfile)
	# load file
	tmp_npy = np.load(filepath)
	if 'merge_array' not in locals():
		ntimes = tmp_npy.shape[0] * n_files
		merge_array = np.zeros((ntimes,tmp_npy.shape[1]))
	else:
		print(start,len(tmp_npy))
		merge_array[start:start+len(tmp_npy)] = copy.copy(tmp_npy)
		start = start+len(tmp_npy)
merge_array_reduced = np.delete(merge_array,np.arange(start,ntimes),axis=0)
np.savetxt(args.output,merge_array_reduced)
np.save(args.output,merge_array_reduced)
print(" total merged npy file = {}".format(merge_array_reduced.shape))