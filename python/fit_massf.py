#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 12/04/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='fitting density profile with tanh and erf function')
## args
parser.add_argument('-i', '--input', default='traj.massf.align.avg', nargs='?', 
	help='mass fraction profile (npy file format, exclude .npy)')
parser.add_argument('-g', '--guess', default='CENTER', nargs='?',
	help='initial guess in center value or highest values (CENTER/any)')
parser.add_argument('-symm', '--symmetry', default='YES', nargs='?',
	help='Use symmetry or no symmetry of coexistent mole fractions (YES/any)')
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
massfrac_1d = np.load(args.input)
massfrac_1d = np.transpose(massfrac_1d)
massfrac_1d_avg = massfrac_1d[0]
massfrac_1d_std = massfrac_1d[1]
nbin = len(massfrac_1d_avg)

## fitting functional form
# wr: mole fraction in A-rich phase
# wp: mole fraction in A-poor phase
# b: center of A-rich phase
# 2c: half-width of A-rich phase
# 2lamda: half-width of interface
def tanh_symm(x, wr, b, c, lamda):
	return 1.0-wr+0.50*(2.0*wr-1.0)*(np.tanh((x-b+c)/lamda)-np.tanh((x-b-c)/lamda))
def erf_symm(x, wr, b, c, lamda):
	return 1.0-wr+0.50*(2.0*wr-1.0)*(erf((x-b+c)/lamda)-erf((x-b-c)/lamda))
def tanh_nosymm(x, wr, wp, b, c, lamda):
	return wp+0.50*(wr-wp)*(2.0*wr-1.0)*(np.tanh((x-b+c)/lamda)-np.tanh((x-b-c)/lamda))
def erf_nosymm(x, wr, wp, b, c, lamda):
	return wp+0.50*(wr-wp)*(2.0*wr-1.0)*(erf((x-b+c)/lamda)-erf((x-b-c)/lamda))

## initial guess
if 'CENTER' in args.guess:
	b = int(nbin/2 - 1)
	wr = massfrac_1d_avg[b]
	wp = massfrac_1d_avg[0]
	print("center wr (avg,std) = {} +- {}".format(wr,massfrac_1d_std[b]))
	print("center wp (avg,std) = {} +- {}".format(wp,massfrac_1d_std[0]))
else:
	b = np.argmax(massfrac_1d_avg)
	wr = np.max(massfrac_1d_avg)
	wp = np.min(massfrac_1d_avg)
c = int(nbin/4)
lamda = int(nbin/10)

## curve fit
x_data = np.linspace(1, nbin, num=nbin, endpoint=True)
if 'YES' in args.symmetry:
	tanh_opt, tanh_cov = curve_fit(tanh_symm,x_data,massfrac_1d_avg,p0=[wr,b,c,lamda],sigma=massfrac_1d_std,bounds=(0,[1., nbin, nbin/2., nbin/2.]))
	erf_opt, erf_cov = curve_fit(erf_symm,x_data,massfrac_1d_avg,p0=[wr,b,c,lamda],sigma=massfrac_1d_std,bounds=(0,[1., nbin, nbin/2., nbin/2.]))
else:
	tanh_opt, tanh_cov = curve_fit(tanh_nosymm,x_data,massfrac_1d_avg,p0=[wr,wp,b,c,lamda],sigma=massfrac_1d_std,bounds=(0,[1., 1., nbin, nbin/2., nbin/2.]))
	erf_opt, erf_cov = curve_fit(tanh_nosymm,x_data,massfrac_1d_avg,p0=[wr,wp,b,c,lamda],sigma=massfrac_1d_std,bounds=(0,[1., 1., nbin, nbin/2., nbin/2.]))

## plotting
if 'YES' in args.show:
	plt.plot(x_data, massfrac_1d_avg, 'b-', label='data')
	if 'YES' in args.symmetry:
		plt.plot(x_data, tanh_symm(x_data,*tanh_opt), 'r--',label='fit:tanh_symm')
		plt.plot(x_data, erf_symm(x_data,*erf_opt), 'g--',label='fit:erf_symm')
	else:
		plt.plot(x_data, tanh_nosymm(x_data,*tanh_opt), 'r--',label='fit:tanh_nosymm')
		plt.plot(x_data, erf_nosymm(x_data,*erf_opt), 'g--',label='fit:erf_nosymm')
	plt.legend()
	#plt.show()
	plt.savefig(args.output_png)

## display all information
if 'YES' in args.symmetry:
	print("tanh wr = {} +- {}".format(tanh_opt[0],tanh_cov[0][0]))
	print("tanh b = {} +- {}".format(tanh_opt[1],tanh_cov[1][1]))
	print("tanh c = {} +- {}".format(tanh_opt[2],tanh_cov[2][2]))
	print("tanh lamda = {} +- {}".format(tanh_opt[3],tanh_cov[3][3]))
	print("erf wr = {} +- {}".format(erf_opt[0],erf_cov[0][0]))
	print("erf b = {} +- {}".format(erf_opt[1],erf_cov[1][1]))
	print("erf c = {} +- {}".format(erf_opt[2],erf_cov[2][2]))
	print("erf lamda = {} +- {}".format(erf_opt[3],erf_cov[3][3]))
else:
	print("tanh wr = {} +- {}".format(tanh_opt[0],tanh_cov[0][0]))
	print("tanh wp = {} +- {}".format(tanh_opt[1],tanh_cov[1][1]))
	print("tanh b = {} +- {}".format(tanh_opt[2],tanh_cov[2][2]))
	print("tanh c = {} +- {}".format(tanh_opt[3],tanh_cov[3][3]))
	print("tanh lamda = {} +- {}".format(tanh_opt[4],tanh_cov[4][4]))
	print("erf wr = {} +- {}".format(erf_opt[0],erf_cov[0][0]))
	print("erf wp = {} +- {}".format(erf_opt[1],erf_cov[1][1]))
	print("erf b = {} +- {}".format(erf_opt[2],erf_cov[2][2]))
	print("erf c = {} +- {}".format(erf_opt[3],erf_cov[3][3]))
	print("erf lamda = {} +- {}".format(erf_opt[4],erf_cov[4][4]))

## timer
hjung.time.end_print(start_proc, start_prof)