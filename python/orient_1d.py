#!/usr/bin/env python3
# ver 0.1 - make codes on 3/29/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculation orientational parameters using com')
## args
parser.add_argument('-icell', '--icell', default='unit_cell.npy', nargs='?', 
	help='input unit cell dimension file')
parser.add_argument('-icom', '--icom', default='pol.com.npy', nargs='?', 
	help='input COM file')
parser.add_argument('-ivec', '--ivec', default='pol.ree.vec.npy', nargs='?', 
	help='input vector file')
parser.add_argument('-nmol', '--nmol', nargs='?', type=int,
	help='# molecules')
parser.add_argument('-axis', '--axis', default=2, nargs='?', type=int,
	help='axis for distribution')
parser.add_argument('-nbin', '--nbin', nargs='?', type=int,
	help='#bins for distribution on a given axis (should be matched with nbins when convolution alignment did)')
parser.add_argument('-o', '--output', default='pol', nargs='?', 
	help='output prefix filename for oriental paramter files (.orient)')

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

args.output = args.output + '.orient'

## read files
unitcell = hjung.io.read_simple(args.icell,0,-1)
n_frames = len(unitcell)
com = hjung.io.read_simple(args.icom,0,-1)
if n_frames != int(len(com)/args.nmol/3):
	raise ValueError("may be wrong n_frames of com data")
else:
	com = com.reshape(n_frames,args.nmol,3)
ree_vec = hjung.io.read_simple(args.ivec,0,-1)
if n_frames != int(len(ree_vec)/args.nmol/3):
	raise ValueError("may be wrong n_frames of ree_vec data")
else:
	ree_vec = ree_vec.reshape(-1,3)

# calc. com histograms
unit_cells_1d   = unitcell[:,args.axis]
com_1d = com[:,:,args.axis]
com_hist_1d_t, bin_1d_t = hjung.analyze.histo_t_1d_nbin(com_1d, unit_cells_1d, args.nbin) 

# calc. orient. 
# see following reference:
#  Structural and thermodynamic properties of interfaces between coexisting phases in polymer blends: a Monte Carlo simulation
#  Marcus Müller, Kurt Binder and Wilfried Oed,  J. Chem. Soc., Faraday Trans., 1995, 91, 2369
#  DOI: 10.1039/FT9959102369
# the variable "orient" is the same as eq.(28) in the paper above.
orient_abs = LA.norm(ree_vec,axis=1)**2 
ree_p = ree_vec[:,args.axis]**2
ree_l1 = ree_vec[:,(args.axis-1)%3]**2
ree_l2 = ree_vec[:,(args.axis-2)%3]**2
orient = (ree_p - 0.5*(ree_l1+ree_l2))/orient_abs
orient = orient.reshape(n_frames,args.nmol)

# calc. average orient. in bins
orient_histo_1d_t = hjung.analyze.histo_xy_t_1d_wbin(com_1d, orient, bin_1d_t)

# save raw rg data file
np.savetxt(args.output, orient_histo_1d_t, 
	header='orientational parameters' , fmt='%f', comments='# ')
np.save(args.output, orient_histo_1d_t)
print(" saved orient_histo files")

## timer
hjung.time.end_print(start_proc, start_prof)