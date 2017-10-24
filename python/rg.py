#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 10/23/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='make index file for radius gyration of selected polymers with fixed degree of polymerization, N_pol')
## args
parser.add_argument('-b', '--begin',  default=1, nargs='?', type=int,
        help='index of beginning atom in conf.gro')
parser.add_argument('-e', '--end', default=-1, nargs='?', type=int,
        help='index of end atom in conf.gro.')
parser.add_argument('-n', '--npol', default=1, nargs='?', type=int,
	help='degree of polymerization of each polymer chains')
parser.add_argument('-o', '--output', default='index.ndx', nargs='?', 
	help='output prefix for radius gyration trajectory')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# check args
print(" input arguments: {0}".format(args))

## import modules
import hjung
from hjung import *
import numpy as np

## timer
start_proc, start_prof = hjung.time.init()

file = open(args.output, 'w')
print("="*30)
init = args.begin
imol = 0
tmol = (args.end - init)/args.npol
while imol < tmol:
	file.write("[ i{} ]\n".format(imol))
	line = list(range(init, init+args.npol))
	file.write(" ".join(map(str,line)))
	file.write("\n")
	file.write("\n")
	imol = imol + 1
	init = init + args.npol
file.close()

## timer
hjung.time.end_print(start_proc, start_prof)