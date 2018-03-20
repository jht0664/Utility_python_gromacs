#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 03/19/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='fitting local density profile with gaussian function')
## args
parser.add_argument('-i', '--input', default='traj.tmass.align.avg', nargs='?', 
	help='local density profile (npy file format, exclude .npy)')
parser.add_argument('-g', '--guess', default='QUARTER', nargs='?',
	help='initial guess in the value in quarter position or least values (QUARTER/any)')
parser.add_argument('-symm', '--symmetry', default='YES', nargs='?',
	help='Use symmetry or no symmetry of coexistent local density (YES/any)')
parser.add_argument('-show', '--show', default='YES', nargs='?', 
	help='Save plotting (YES/any)')
parser.add_argument('-o', '--output', default='.fit', nargs='?', 
	help='output surfix for fitting result')
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
from scipy.special import erf
from scipy.optimize import curve_fit
import math
import matplotlib
matplotlib.use('Agg') # avoid to show figures when running bash shell script
import matplotlib.pyplot as plt

# default for args
args.input = args.input + '.npy'
args.output = args.input + args.output
args.output_png = args.output + '.png'

## timer
start_proc, start_prof = hjung.time.init()

## load data files
tmass_1d = np.load(args.input)
#tmass_1d = np.loadtxt(args.input)
tmass_1d = np.transpose(tmass_1d)
tmass_1d_avg = tmass_1d[0]
tmass_1d_std = tmass_1d[1]
curve_fit_std_off = False
if len(np.nonzero(tmass_1d_std)) != len(tmass_1d_std):
	print("local density std elements have zeros. Turned off curve_fit using std.")
	curve_fit_std_off = True
nbin = len(tmass_1d_avg)

## fitting normal distribution form (gaussian)
def gaus_symm(x, pref, sigma, itf_pos1, itf_pos2, bulk_tmass):
	return -pref*np.exp(-(x-itf_pos1)**2.0/(2.0*sigma**2.0)) - pref*np.exp(-(x-itf_pos2)**2.0/(2.0*sigma**2.0)) + bulk_tmass

## initial guess
if 'QUARTER' in args.guess:
	print("initial guess algorithm is active")
	itf_pos1 = int(nbin/4 - 1)
	itf_pos2 = int(itf_pos1 + nbin/2.)
	center = int(nbin/2 - 1)
	pref = tmass_1d_avg[itf_pos1]
	bulk_tmass = tmass_1d_avg[center]
else:
	itf_pos1 = np.argmin(tmass_1d_avg)
	itf_pos2 = int(itf_pos1 + nbin/2.)
	pref = np.min(tmass_1d_avg)
	bulk_tmass = np.max(tmass_1d_avg)
gaus_var = nbin/4.0
#gaus_var = (2.0*math.pi*pref**2.0)**(-0.50)
print(" initial guess = {} pref, {} gaus_var, {} itf_pos1, {} itf_pos2, {} bulk_tmass".format(pref,gaus_var,itf_pos1,itf_pos2,bulk_tmass))

## curve fit
x_data = np.linspace(1, nbin, num=nbin, endpoint=True)
if 'YES' in args.symmetry:
	gaus_opt, gaus_cov = curve_fit(gaus_symm,x_data,tmass_1d_avg,
		p0=[pref,gaus_var,itf_pos1,itf_pos2,bulk_tmass],sigma=tmass_1d_std,
		bounds=([0,0,0,nbin/2.,0],[100, nbin/4., nbin/2., nbin, 100])) # maybe option bounds makes error to fit
else:
	raise ValueError("not supported yet for asymmetric profile")
	#if curve_fit_std_off:
		#tanh_opt, tanh_cov = curve_fit(tanh_nosymm,x_data,tmass_1d_avg,p0=[wr,wp,b,c,lamda],bounds=(0,[1., 1., nbin, nbin/2., nbin/2.]))
		#erf_opt, erf_cov = curve_fit(tanh_nosymm,x_data,tmass_1d_avg,p0=[wr,wp,b,c,lamda],bounds=(0,[1., 1., nbin, nbin/2., nbin/2.]))
	#else:
		#tanh_opt, tanh_cov = curve_fit(tanh_nosymm,x_data,tmass_1d_avg,p0=[wr,wp,b,c,lamda],sigma=tmass_1d_std,bounds=(0,[1., 1., nbin, nbin/2., nbin/2.]))
		#erf_opt, erf_cov = curve_fit(tanh_nosymm,x_data,tmass_1d_avg,p0=[wr,wp,b,c,lamda],sigma=tmass_1d_std,bounds=(0,[1., 1., nbin, nbin/2., nbin/2.]))

## plotting
if 'YES' in args.show:
	plt.plot(x_data, tmass_1d_avg, 'b-', label='data')
	if 'YES' in args.symmetry:
		plt.plot(x_data, gaus_symm(x_data,*gaus_opt), 'r--',label='fit:gaussian_symm')
	#else:
		#plt.plot(x_data, tanh_nosymm(x_data,*tanh_opt), 'r--',label='fit:tanh_nosymm')
		#plt.plot(x_data, erf_nosymm(x_data,*erf_opt), 'g--',label='fit:erf_nosymm')
	plt.legend()
	#plt.show()
	plt.savefig(args.output_png)

## display all information
if 'YES' in args.symmetry:
	print("norm. factor (peak value) = {} +- {}".format(gaus_opt[0],gaus_cov[0][0]))
	#print("norm. var = {} +- {}".format(gaus_opt[1],gaus_cov[1][1]))
	print("interface pos1 = {} +- {}".format(gaus_opt[2],gaus_cov[2][2]))
	print("interface pos2 = {} +- {}".format(gaus_opt[3],gaus_cov[3][3]))
	print("bulk tmass = {} +- {}".format(gaus_opt[4],gaus_cov[4][4]))
#else:
	#print("tanh wr = {} +- {}".format(tanh_opt[0],tanh_cov[0][0]))
	#print("tanh wp = {} +- {}".format(tanh_opt[1],tanh_cov[1][1]))
	#print("tanh b = {} +- {}".format(tanh_opt[2],tanh_cov[2][2]))
	#print("tanh c = {} +- {}".format(tanh_opt[3],tanh_cov[3][3]))
	#print("tanh lamda = {} +- {}".format(tanh_opt[4],tanh_cov[4][4]))
	#print("erf wr = {} +- {}".format(erf_opt[0],erf_cov[0][0]))
	#print("erf wp = {} +- {}".format(erf_opt[1],erf_cov[1][1]))
	#print("erf b = {} +- {}".format(erf_opt[2],erf_cov[2][2]))
	#print("erf c = {} +- {}".format(erf_opt[3],erf_cov[3][3]))
	#print("erf lamda = {} +- {}".format(erf_opt[4],erf_cov[4][4]))

## timer
hjung.time.end_print(start_proc, start_prof)