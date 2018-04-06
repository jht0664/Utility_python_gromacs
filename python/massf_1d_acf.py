#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 04/06/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='autocorrelation functions along time for a given position')
## args
parser.add_argument('-i', '--input', default='peo.massf.align.npy', nargs='?', 
	help='input massf or 2d file')
parser.add_argument('-o', '--output', default='.time_acf', nargs='?', 
	help='output prefix 2d acf file along time')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import sys
sys.path.append('/home/htjung/Utility/python/')
import hjung
from hjung import *
import numpy as np
from scipy import ndimage

# default for args
args.output = args.output + '.time_acf'

## timer
start_proc, start_prof = hjung.time.init()

# load data
data_1d_t = np.load(args.input)
n_frames = len(data_1d_t)
data_t_1d = np.transpose(data_1d_t)
nbins = len(data_t_1d)

del_data = data_t_1d - np.mean(data_t_1d,axis=1).reshape((nbins,-1))
norm = np.sum(del_data**2,axis=1)
acf_data = np.zeros((nbins,n_frames))
for i in range(nbins):
	acf_data[i] = np.correlate(del_data[i],del_data[i],mode='same')/norm[i]
acf_1d_t = np.delete(np.transpose(acf_data),np.arange(int(n_frames/2.)),axis=0)

## save number histogram trajectory
np.savetxt(args.output, acf_1d_t, 
	header='{} autocorrelations for 2d profile along time at columns (space)'.format(nbins), fmt='%f', comments='# ')
np.save(args.output, acf_1d_t)

## timer
hjung.time.end_print(start_proc, start_prof)