#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 9/08/2016

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Arimethic mean and standard deviation')

# args
parser.add_argument('-i', '--input', nargs='?', help='input file')
#parser.add_argument('-o', '--output', nargs='?', help='output file')
parser.add_argument('-b', '--begin', type=int, default=0, nargs='?', 
	help='# of lines to ignore for statistics')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()

# default
args.begin = args.begin if args.begin is not None else 0

# Check arguments for log
print("input filename   = ", args.input)
print("# skipped lines    = ", args.begin)
#print("output filename  = ", args.output)

# Read data
try:
	xvg = open(args.input, 'r')
except IOError:
	print("Problem with opening ",args.input)
	exit()
data = []
i = 0
import numpy as np
for line in xvg:
	i += 1
	if i <= args.begin:
		continue
	line = line.strip()
	if not line or line.startswith('#') or line.startswith('@'): 
		# line is blank or comment line
		continue
	data = np.append(data, np.float64(line))
xvg.close()

# statistics
print("======================================")
print("Arimetic mean      = ", np.mean(data))
print("Standard deviation = ", np.std(data))
