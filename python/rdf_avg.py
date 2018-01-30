#!/usr/bin/env python3
# ver 0.1 - copy from rdf_itf.py (v0.1) and modify codes on 1/29/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='average rdf(s)')
## args
parser.add_argument('-mi', '--max_i', default=10, nargs='?', type=int,
	help='max number of folder')
parser.add_argument('-ir', '--ir', default=0, nargs='?', type=int,
	help='numbering of partial rdfs to avg')
parser.add_argument('-folder', '--folder_name', default='nvt', nargs='?', 
	help='surfix foldername (like *nvt)')
parser.add_argument('-file', '--file_name', default='a.massf.domain.rdf', nargs='?',
	help='prefix of rdfs filename')
parser.add_argument('-o', '--output', default='.avg', nargs='?', 
	help='output prefix for avg rdf files')
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

# default for args
oavg = args.file_name + '.' + str(args.ir) + args.output

## timer
start_proc, start_prof = hjung.time.init()

## load rdf files
for i in range(args.max_i):
	rdf_tmp = np.load(str(i+1)+args.folder_name+'/'+args.file_name+'.'+str(args.ir)+'.npy')
	bin_edges = rdf_tmp[:,0]
	hist_values = rdf_tmp[:,1]
	if i == 0:
		rdf_avg = np.zeros(len(hist_values))
		rdf_std = np.zeros(len(hist_values))
	rdf_avg = rdf_avg + hist_values
	#print(rdf_avg)
	rdf_std = rdf_std + (hist_values**2.0)
	#print(rdf_std)

rdf_avg = rdf_avg/float(args.max_i)
rdf_std = rdf_std/float(args.max_i)
rdf_std = rdf_std - np.square(rdf_avg)
rdf_data = np.column_stack((bin_edges,rdf_avg, np.sqrt(rdf_std)))

np.savetxt(oavg, rdf_data, 
	header='averaged radial distribution functions', fmt='%e', comments='# ')
np.save(oavg, rdf_data)

## timer
hjung.time.end_print(start_proc, start_prof)