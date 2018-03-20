#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 03/19/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Calculate distribution of surface tension from normal and traverse pressures from test area methods')
## args
parser.add_argument('-ipn', '--ipn', default='surf_dist_z_z.out', nargs='?', 
	help='normal pressure results by test area methods')
parser.add_argument('-ipt', '--ipt', default='surf_dist_xy_z.out', nargs='?', 
	help='traverse pressure results by test area methods')
parser.add_argument('-o', '--surfix', default='.avg', nargs='?', 
	help='surfix of output avg file.')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np

## timer
start_proc, start_prof = hjung.time.init()

## output filename
args.opn = args.ipn + args.surfix
args.opt = args.ipt + args.surfix

## load data files
raw_pn = np.loadtxt(args.ipn,comments='#',unpack=True)
raw_pt = np.loadtxt(args.ipt,comments='#',unpack=True)
if raw_pn.size != raw_pt.size:
	raise ValueError("the size of two data files are different.")
nwindows = len(raw_pn)
nframes = len(raw_pn[0])

## calculation
raw_pn_avg = np.average(raw_pn,axis=1)
raw_pn_std = np.std(raw_pn,axis=1)
raw_pt_avg = np.average(raw_pt,axis=1)
raw_pt_std = np.std(raw_pt,axis=1)

# for graphing
raw_pn_dataset = np.column_stack((raw_pn_avg,raw_pn_std))
raw_pt_dataset = np.column_stack((raw_pt_avg,raw_pt_std))

## fit
from scipy import stats

x = np.arange(1,nwindows+1)
pn, pn_cov = np.polyfit(x,raw_pn_avg,1,cov=True) # pn[0]: slope, pn[1]: y-intercept
a, b, pn_r, c, d = stats.linregress(x,raw_pn_avg)
print("P_N (R2) = {} ".format(pn_r**2)) 
pt, pt_cov = np.polyfit(x,raw_pt_avg,1,cov=True)
a, b, pt_r, c, d = stats.linregress(x,raw_pt_avg)
print("P_T (R2) = {} ".format(pt_r**2))

print("P_N = {} +- {}".format(pn[1], np.sqrt(pn_cov[1][1]))) 
print("P_T = {} +- {}".format(pt[1], np.sqrt(pt_cov[1][1]))) 
gamma_cov = pn_cov[1][1] + pt_cov[1][1] # see details about adding/substrating two independnt numbers 
# https://stats.stackexchange.com/questions/112351/standard-deviation-after-subtracting-one-mean-from-another
print("gamma/Lz = {} +- {}".format((pn[1] - pt[1])/2.0, np.sqrt(gamma_cov)/2.0)) 

## write
np.savetxt(args.opn, raw_pn_dataset, 
	header='avg std for normal pressure with %d windows' \
	%(nwindows), fmt='%f', comments='# ')
np.savetxt(args.opt, raw_pt_dataset, 
	header='avg std for parallel pressure with %d windows' \
	%(nwindows), fmt='%f', comments='# ')

## timer
hjung.time.end_print(start_proc, start_prof)