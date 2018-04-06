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

start = 0
with open(args.output,'wb') as f_handle:
	npyfiles = glob.glob("*.npy")
	npyfiles.sort()
	n_files = len(npyfiles)
	for npyfile in npyfiles:
		# get path
		filepath = os.path.join('./', npyfile)
		# load file
		tmp_npy = np.load(filepath)
		if 'merge_array' not in locals():
			ntimes = tmp_npy.shape[0] * n_files
			merge_array = np.zeros((ntimes,tmp_npy.shape[1]))
		else:
			merge_array[start:start+len(tmp_npy)] = copy.copy(tmp_npy)
			start = start+len(tmp_npy)
	np.save(f_handle,merge_array)
print(merge_array.shape)