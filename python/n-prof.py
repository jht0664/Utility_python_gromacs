#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 11/15/2016

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='1D Density profile')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr structure file')
parser.add_argument('-select', '--select', nargs='?', 
	help='a file with a command-line for select_atoms in MDAnalysis')
parser.add_argument('-bin', '--bin', default=2.0, nargs='?', 
	help='bin size (A), defult 2A')
parser.add_argument('-nbin', '--nbin', nargs='?', 
	help='number of bins, otherwise we use the bin size')
parser.add_argument('-axis', '--axis', default=2, nargs='?', 
	help='which axis for histogram (x axis (0), y axis (1), z axis (2))')
parser.add_argument('--nojump', action='store_true',
	help='re-positioning (x,y,z) within an unit cell (option)')
parser.add_argument('-onum', '--outputnumber', default='traj.nums', nargs='?', 
	help='output file for number trajectory')
parser.add_argument('-obin', '--outputbin', default='traj.bins', nargs='?', 
	help='output file for bin trajectory')
parser.add_argument('-o', '--output', default='traj.nprob', nargs='?', 
	help='output file for number probabilty')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# default for args
args.input = args.input if args.input is not None else 'traj.trr'
args.structure = args.structure if args.structure is not None else 'topol.tpr'
args.output = args.output if args.output is not None else 'traj.nprob'
args.outputnumber = args.outputnumber if args.outputnumber is not None else 'traj.nums'
args.outputbin = args.outputbin if args.outputbin is not None else 'traj.bins'
args.bin = args.bin if args.bin is not None else 2.0
args.axis = args.axis if args.axis is not None else 2
args.bin = float(args.bin)
args.axis = int(args.axis)
args.nbin = int(args.nbin)

## Check arguments for log
print("===============================")
print("input filename   = ", args.input)
print("str filename     = ", args.structure)
if args.input is not None:
	print("select filename  = ", args.select)
print("No Jump if yes   = ", args.nojump)
print("bin size (A)     = ", args.bin)
print("number of bins   = ", args.nbin)
print("axis [0:2]       = ", args.axis)
print("output number trajectory filename = ", args.outputnumber)
print("output bin trajectory filename = ", args.outputbin)
print("output number probability filename  = ", args.output)

## check vaulable setting
if args.bin < 1.0:
	raise ValueError("Too small bin size is usually not valid for histogram")
if args.axis < 0 or args.axis > 2:
	raise ValueError("wrong input of axis for histogram")
if args.nbin is not None:
	print("BE AWARE: You use number of bins, instead of bin size")
	use_nbin = True
else:
	use_nbin = False

## timer
import time
start_clock = time.clock() # process time
start_wall = time.time() # wall time

# Read 3D-coordinates from .trr trajectory files 
# input: tpr_filename, the filename of .tpr file in Gromacs
#		 trr_filename, the filename of .trr file in Gromacs
#		 select_atoms_filename, the filename including a command-line
#			for keeping trajectory of the selected atoms
# output: coordinates, xyz position of atoms, (#frame x #atoms x 3)
#		  unit_cells, box dimension (#frame x (x,y,z)) 
def read_coord_trr_3d(tpr_filename, trr_filename, select_atoms_filename):
	# read a line of select command-line in MDAnalysis
	if select_atoms_filename is not None:
		try:
			open_file = open(select_atoms_filename, 'r')
		except IOError:
			raise IOError("Problem with opening ",select_atoms_filename)
		print("reading the first line in %s file" %select_atoms_filename)
		select_command = open_file.readline().strip()
		print("select: %s" %select_command)
		open_file.close()
	# import
	import MDAnalysis
	import numpy as np
	# set up for reading trajectory
	u = MDAnalysis.Universe(tpr_filename,trr_filename)
	if select_atoms_filename is None:
		n_atoms = len(u.atoms)
	else:
		atoms = u.select_atoms(select_command).indices
		n_atoms = len(atoms)
	print("# atoms  = ", n_atoms)
	n_frames = len(u.trajectory)
	print("# frames = ", n_frames)
	# initailize variables
	coordinates = np.zeros((n_frames, n_atoms, 3))
	#velocities = np.zeros(n_frames, n_atoms, 3)
	#forces = np.zeros(n_frames, n_atoms, 3)
	unit_cells = np.zeros((n_frames, 6))
	# read trajectory
	i_frame = 0
	for ts in u.trajectory:
		try:
			if select_atoms_filename is None:
				coordinates[i_frame, :, :] = ts._pos
			else: 
				tmp = np.array(ts._pos)
				coordinates[i_frame, :, :] = tmp[atoms] 
				# read positions in 10^-10 m (A)
			#np.append(velocities, ts._velocities)
			#np.append(forces, ts._forces)
			unit_cells[i_frame, :] = ts._unitcell
		except IndexError:
			raise ValueError("There are more coordinates to be read than indicated in the header.")
		i_frame += 1
	# check consistency; final i_frame should be the same as # frames
	if i_frame != n_frames:
		raise ValueError("# of frames to read %d does not agree with the length of trajectory file %d" \
	                             % (i_frame, n_frames))
	# box info
	if all(unit_cells[0,:] == unit_cells[1,:]):
		print("The system may be in NVT ensemble")
	else:
		if unit_cells[0][0] == unit_cells[1][0] and unit_cells[0][1] == unit_cells[1][1]:
			print("may be in NPAT ensemble")
		else:
			print("may be in NPT ensemble")
	return coordinates, unit_cells
## read a topology and a trajectory using module MDAnalysis with selection
print("===============================")
coordinates, unit_cells = read_coord_trr_3d(args.structure, args.input, args.select)
print("Done: reading trajectory and topology file")

# move atoms inside simulation box (no jump option)
# input: xyz_t is coordinate sets of atoms along time (t1, t2, ...)
#			[[[x, y, z], [x, y, z], ...], [[x, y, z], [x, y, z], ...], ...]
# 		 box_t is box dimension sets along time (t1, t2, ...)
#			[[[box_x, box_y, box_z]], [[box_x, box_y, box_z]], ...]
# output: a new confined coordinate inside box along time
def pbc_nojump_t(xyz_t, box_t): 
	import numpy as np
	# convert double float
	xyz_t = np.array(xyz_t, dtype=np.float64)
	box_t = np.array(box_t, dtype=np.float64)
	# start for loop
	print("no jump option activated")
	i_frame = 0;
	for time in xyz_t:
		# check if the unit cell is rectangular parallelepiped
		if box_t[i_frame][3] != box_t[i_frame][4] or box_t[i_frame][4] != box_t[i_frame][5] or box_t[i_frame][3] != box_t[i_frame][5]:
			raise TypeError("Not support except rectangular parallelepiped, ", box_t[i_frame][3:5])
		# wrap coodinates within unit cell
		for atom_xyz in time:
			for xyz in atom_xyz:
				xyz[0] = xyz[0] - box_t[i_frame][0]*np.floor(xyz[0]/box_t[i_frame][0])# x position
				xyz[1] = xyz[1] - box_t[i_frame][1]*np.floor(xyz[1]/box_t[i_frame][1])# y position 
				xyz[2] = xyz[2] - box_t[i_frame][2]*np.floor(xyz[2]/box_t[i_frame][2])# z position 
		i_frame += 1
	return xyz_t
## wrap positions within unit cell
if args.nojump:
	print("===============================")
	coordinates = pbc_nojump_t(coordinates, unit_cells)
	print("Done: nojumping coordinates")

## reduce 3d-coordinates to 1d-coordinates
coordinates_1d = coordinates[:,:,args.axis]
unit_cells_1d = unit_cells[:,args.axis]

# make 1D histograms every single frame using 1d-coordinates
# input: x_t is 1d-coordinate of atoms along time (t1, t2, ...)
#			[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
# 		 box_t is box dimension sets along time (t1, t2, ...)
#			[box_x(t0), box_x(t1), ...]
#        bin_size is which axis you want to do histogram
#		 nbin is the number of bins you want 
#		 use_nbin is true if you want to use nbin instead of bin_size. otherwise, false
# output: histo_t is a new 1d-number profile trajectory
#		  bin_t is a bin position array for the 1d histogram 
def histo_t_1d(x_t, box_x_t, bin_size, nbin, use_nbin):
	import numpy as np
	# check n_frames of x_t and box_x_t
	if len(x_t) != len(box_x_t):
		raise ValueError("# of time frame is not the same for input arrarys")
	# set number of bins (const) using input and initial box dimension
	if use_nbin:
		print("bin size = %d", box_x_t[0]/float(nbin))
		n_bins = nbin
	else:
		n_bins = int(np.around(box_x_t[0]/bin_size))
	# initailize variables
	n_frames = len(x_t)
	histo_t = np.zeros((n_frames, n_bins))
	bin_t = np.zeros((n_frames, n_bins+1))
	# make histogram trajectory
	i_frame = 0
	for x in x_t:
		histo_t[i_frame], bin_t = np.histogram(x, bins=n_bins, range=(0,box_x_t[i_frame]))
		i_frame += 1
	return histo_t, bin_t

## number histograms for each frame 
print("===============================")
number_t_1d, bin_t_1d = histo_t_1d(coordinates_1d, unit_cells_1d, args.bin, args.nbin, use_nbin) 
print("Done: making number trajectory with respect to bins")

## save number histograms
print("===============================")
import numpy as np
np.savetxt(args.outputnumber, number_t_1d, fmt='%d')
np.savetxt(args.outputbin, bin_t_1d)
print("Finished saving files of number and bin trajectory")

# merge all arraylists to one arraylist in order to remove time info.
# input: x_t is array along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
# output: [x1(t0), x2(t0), x3(t0), ..., x1(t1), x2(t1), x3(t1), ...] 
def merge_t_to_1d(x_t):
	# internally merge arraylist in a list
	import itertools
	# merge 
	return [i for i in itertools.chain.from_iterable(x_t)]

## make an array for probability of numbers and histogram, then save the probability in a file
number_1d = merge_t_to_1d(number_t_1d)

# make 1d density histrogram of numbers
# input: x is 1d-array
# output: [[a collection set of density of each bins],[a collection set of bin]]]
def histo_num_dens_1d(x):
	import numpy as np
	n_bins = int(max(x))
	return np.histogram(x, bins=n_bins, density=True)

## make 1d density probability of number
print("===============================")
prob_num_1d, prob_num_bin_1d = histo_num_dens_1d(number_1d)
print("Done: making a time- and bin-averaged probability of numbers")

## re-assign bin value for writing a file
print("===============================")
for i in range(1,len(prob_num_bin_1d)-1):
	prob_num_bin_1d[i-1] = 0.50*(prob_num_bin_1d[i-1]+prob_num_bin_1d[i])
prob_num_bin_1d = prob_num_bin_1d[0:len(prob_num_bin_1d)-1]
## save file
np.savetxt(args.output, np.column_stack((prob_num_bin_1d, prob_num_1d)), fmt="%7.3f %0.8f")
print("Finished saving files for probability of number")

## translate to get the lowest RMSE

## save the density profile trajectory

## block average

## save the block average

print("===============================")
print(time.clock() - start_clock, "seconds process time")
print(time.time() - start_wall, "seconds wall time")