#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 2/10/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Generate atom (index) pairs of PEO/PEG for Gromacs ')
## args
parser.add_argument('-deg', '--degree', nargs='?', 
	help='Degree of polymerization, n, of R1-[CH2-O-CH2]n-R1')
parser.add_argument('-r1', '--r1', nargs='?', 
	help='Beginning group, R1, CH2OH for PEG')
parser.add_argument('-r2', '--r2', nargs='?', 
	help='End group, R2, CH3OCH3 for PEG/PEO')
parser.add_argument('-o', '--output', default='bond.pair', nargs='?', 
	help='output file for bond pairs for topology file of PEG/PEO')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# default for args
args.degree = int(args.degree) 
args.output = args.output if args.output is not None else 'bond.pair'

## Check arguments for log
print("===============================")
print("begining group  = ", args.r1)
print("end      group  = ", args.r2)
print("degree of polymerization = ", args.degree)
print("PEG/PEO is R1-[CH2-O-CH2]n-R1")
print("output filename = ", args.output)

#Estimate total natoms
natoms = args.degree*7
if args.r1 == 'CH2OH':
	natoms = natoms + 5
if args.r2 == 'CH2OCH3':
	natoms = natoms + 8
print("total estimated number of atoms = ", natoms)

with open(args.output, 'w') as output_file:
	
	output_file.write("%d \n" %natoms)
	# Bond
	index = 0
	output_file.write("[ bonds ] \n")
	output_file.write("; ai  aj  funct \n")
	print("Bond: beginning group, %s" %(args.r1))
	if args.r1 == 'CH2OH':
		output_file.write(" %d   %d   1\n" %(index+1,index+2))
		output_file.write(" %d   %d   1\n" %(index+2,index+3))
		output_file.write(" %d   %d   1\n" %(index+3,index+4))
		output_file.write(" %d   %d   1\n" %(index+3,index+5))
		output_file.write(" %d   %d   1\n" %(index+3,index+6))
		index = index + 6
	else:
		raise ValueError("Not supported yet for beginning group, %s" %args.r1)
	
	print("Bond: [-CH2-O-CH2-]n chain")
	if args.degree < 3:
		raise ValueError("Not supported yet for short chain (%d-mer)" %args.degree)
	for deg in range(args.degree):
		output_file.write(" %d   %d   1\n" %(index,index+1))
		output_file.write(" %d   %d   1\n" %(index,index+2))
		output_file.write(" %d   %d   1\n" %(index,index+3))
		output_file.write(" %d   %d   1\n" %(index+3,index+4))
		output_file.write(" %d   %d   1\n" %(index+4,index+5))
		output_file.write(" %d   %d   1\n" %(index+4,index+6))
		output_file.write(" %d   %d   1\n" %(index+4,index+7))
		index = index + 7

	print("Bond: a end group, %s" %(args.r2))
	if args.r2 == 'CH2OCH3':
		output_file.write(" %d   %d   1\n" %(index,index+1))
		output_file.write(" %d   %d   1\n" %(index,index+2))
		output_file.write(" %d   %d   1\n" %(index,index+3))
		output_file.write(" %d   %d   1\n" %(index+3,index+4))
		output_file.write(" %d   %d   1\n" %(index+4,index+5))
		output_file.write(" %d   %d   1\n" %(index+4,index+6))
		output_file.write(" %d   %d   1\n" %(index+4,index+7))
		index = index + 7
	else:
		raise ValueError("Not supported yet for end group, %s" %args.r2)

	print("Total # atoms = %d" %index)
	if index != natoms:
		raise ValueError("Estimation is wrong, Need to check agina")
	output_file.write("\n")
	output_file.write("\n")

	print("Finish!")







