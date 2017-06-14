#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Calculate instantaneous surface tension from pressure tensor')
# args
parser.add_argument('-box', '--box', type=float, nargs=3, 
	help='x, y, z (nm) of box dimension in NVT')
parser.add_argument('-s', '--structure', nargs='?', 
	help='structure file with box dimension at the last')
parser.add_argument('-i', '--input', default='energy.xvg', nargs='?', 
	help='input file to read pressure tensors [bar], (time, Pxx, Pyy, Pzz)')
parser.add_argument('-o', '--output', default='surfacetension.xvg', nargs='?', 
	help='output file')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# default for args
args.input = args.input if args.input is not None else 'energy.xvg'
args.output = args.output if args.output is not None else 'surfacetension.xvg'

# assign the longest box length
if args.structure is not None and args.box is None:
	import subprocess
	# byte string to string
	lastline = subprocess.check_output(['tail','-1',args.structure]).decode()
	x, y, z = lastline.split()
	x = float(x)
	y = float(y)
	z = float(z)
elif args.structure is None and args.box is not None:
	x = args.box[1]
	y = args.box[2]
	z = args.box[3]
else:
	print("Box dimension is required!")
	exit()
max_box_length = max(x, y, z)


# Check arguments for log
print("the longest box_length = ", max_box_length)
print("input filename         = ", args.input)
print("output filename        = ", args.output)

# Surface tenstion Calculation using pressure tensor
tensor = []
try:
	xvg = open(args.input, 'r')
except IOError:
	print("Problem with opening ",args.input)
	exit()

# surfacetension function
def surfacetension(z_length, Pz, Px, Py):
	return 0.5*z_length*(float(Pz)-0.5*(float(Px)+float(Py)))/10.0

# calculation
for line in xvg:
	line = line.strip()
	if not line or line.startswith('#'): # line is blank or comment line
		continue
	time, xx, yy, zz = line.split()
	# which axis is the longest
	# pressure tensor and box length units are bar and nm in gromacs
	# convert [bar nm] to [mN / m] 
	if max_box_length == x:
		tensor.append(surfacetension(max_box_length,xx,yy,zz))
	elif max_box_length == y:
		tensor.append(surfacetension(max_box_length,yy,xx,zz))
	elif max_box_length == z:
		tensor.append(surfacetension(max_box_length,zz,yy,xx))
	else:
		print("max_box_length does not match box dimension")
		exit()
xvg.close()

# save data
with open(args.output, 'w') as output_file:
	output_file.write("# instantaneous surface tension [mN/m]")
	for item in tensor:
		output_file.write("{}\n".format(item))

# result
print("number of values: ", len(tensor))
