#!/usr/bin/env python3
# ver 0.1 - make codes on 3/29/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculation properties (scalar) distribution using com')
## args
parser.add_argument('-icell', '--icell', default='unit_cell.npy', nargs='?', 
	help='input unit cell dimension file')
parser.add_argument('-icom', '--icom', default='pol.com.npy', nargs='?', 
	help='input COM file')
parser.add_argument('-iprop', '--iprop', default='pol.ree.npy', nargs='?', 
	help='input property (scalar) file')
parser.add_argument('-nmol', '--nmol', nargs='?', type=int,
	help='# molecules')
parser.add_argument('-axis', '--axis', default=2, nargs='?', type=int,
	help='axis for distribution')
parser.add_argument('-nbin', '--nbin', nargs='?', type=int,
	help='#bins for distribution on a given axis (should be matched with nbins when convolution alignment did)')
parser.add_argument('-o', '--output', default='pol.ree', nargs='?', 
	help='output prefix filename for property (.1d)')

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
from numpy import linalg as LA

import MDAnalysis as mda
from MDAnalysis.analysis import distances

from scipy.spatial.distance import euclidean
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot

## timer
start_proc, start_prof = hjung.time.init()

args.output = args.output + '.1d'
ocom = args.icom.replace('.npy','')+'.1d'

## read files
unitcell = hjung.io.read_simple(args.icell,0,-1)
n_frames = len(unitcell)
com = hjung.io.read_simple(args.icom,0,-1)
if n_frames != int(len(com)/args.nmol/3):
	raise ValueError("may be wrong n_frames of com data")
else:
	com = com.reshape(n_frames,args.nmol,3)
prop_scalar = hjung.io.read_simple(args.iprop,0,-1)
if n_frames != len(prop_scalar):
	raise ValueError("may be wrong n_frames of prop_scalar data")
else:
	prop_scalar = prop_scalar.reshape(-1,args.nmol)

# calc. com histograms
unit_cells_1d   = unitcell[:,args.axis]
com_1d = com[:,:,args.axis]
com_hist_1d_t, bin_1d_t = hjung.analyze.histo_t_1d_nbin(com_1d, unit_cells_1d, args.nbin) 
np.savetxt(ocom, com_hist_1d_t, 
	header='com distribution with {} frames and {} bins'.format(n_frames,args.nbin) , fmt='%f', comments='# ')
prop_histo_1d_t = hjung.analyze.histo_xy_t_1d_wbin(com_1d, prop_scalar, bin_1d_t)

# save raw rg data file
np.save(ocom, com_hist_1d_t)
print(" saved com hist files")
np.savetxt(args.output, prop_histo_1d_t, 
	header='property distribution with {} frames and {} bins'.format(n_frames,args.nbin), fmt='%f', comments='# ')
np.save(args.output, prop_histo_1d_t)
print(" saved prop_histo_1d_t files")

## timer
hjung.time.end_print(start_proc, start_prof)