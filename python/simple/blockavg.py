#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 9/08/2016
#         - (Ref.) Appendix D. Statistical Errors 
#       		in the book Understanding Molecular Simulation 
# 				by Daan Frenkel
#		  - (Ref.) Chapter2 in annual Reports in Computational Chemistry, 
#				Vol 5, 2009, doi: 10.1016/S1574-1400(09)00502-7

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Get block average \n'
				'See details in Appendix D. Statistical Errors \n'
				' in the book Understanding Molecular Simulation \n'
				' by Daan Frenkel')
# args
parser.add_argument('-i', '--input', default='block.in', nargs='?', 
	help='input file')
parser.add_argument('-o', '--output', default='block.out', nargs='?', 
	help='output file')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# default for args
args.input = args.input if args.input is not None else 'block.in'
args.output = args.output if args.output is not None else 'block.out'

# Check arguments for log
print("input filename   = ", args.input)
print("output filename  = ", args.output)

# Read data
try:
	xvg = open(args.input, 'r')
except IOError:
	print("Problem with opening ",args.input)
	exit()
data = []
for line in xvg:
	line = line.strip()
	if not line or line.startswith('#') or line.startswith('@'): # line is blank or comment line
		continue
	data.append(float(line))
xvg.close()
print("======================================")
print("# data = ", len(data))

# make a list of block length we will do
import math
list_block_length = []
block_length = 5
while True:
	num_block = math.floor(len(data)/block_length)
	if num_block > 5:
		list_block_length.append(block_length)
		block_length = block_length * 2
		continue
	else:
		break

# computation of block length, # blocks, blocked average,
#  and blocked standard error (= square root of variance = standard deviation)
import statistics # Python3 support this module
import math 
avg = statistics.mean(data)
std = statistics.pstdev(data)
list_bse = []
for block_length in list_block_length:
	list_bavg = []
	for i in range(0,len(data),block_length):
		# mean of each block
		bavg = statistics.mean(data[i:i+block_length-1])
		list_bavg.append(bavg)
	# assume all blocks have the same weight
	bstd = statistics.pstdev(list_bavg)
	bavg_b = statistics.mean(list_bavg)
	bse = bstd/math.sqrt(len(list_bavg))
	list_bse.append([block_length, bse, bavg_b, bstd])

# write list of blocked standard error
import csv
with open(args.output, 'w') as output_file:
	output_file.write("# block length, blocked standard error, blocked average, blocked STD \n")
	output_file.write("# avg: {}, std: {}\n".format(avg, std))
	outputwriter = csv.writer(output_file, delimiter=' ')
	outputwriter.writerows(list_bse)
