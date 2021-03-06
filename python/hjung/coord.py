# move atoms inside simulation box (no jump option)
# input: xyz_t is coordinate sets of atoms along time (t1, t2, ...)
#		[[[x, y, z], [x, y, z], ...], [[x, y, z], [x, y, z], ...], ...]
#	 box_t is box dimension sets along time (t1, t2, ...)
#		[[box_x, box_y, box_z]], [[box_x, box_y, box_z]], ...]
# output: a new confined coordinate inside box along time
# Example:
#if args.nojump:
#        print("===============================")
#        coordinates = pbc_nojump_t(coordinates, unit_cells)
#        print("Done: nojumping coordinates")

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

# wrap atoms within simulation box (no jump option)
# input: x is coordinate sets of atoms
#		[x1, x2, x3, ..]
#	     box is box length
# output: a new wrapped coordinate inside box
def pbc_nojump_1d(x, box_x): 
	import numpy as np
	import copy
	x_out = copy.copy(x)
	#print(" io.pbc_nojump_1d:")
	for ix in range(len(x)):
		x_out[ix] = x[ix] - box_x*np.floor(x[ix]/box_x) # wrap coodinates within unit cell
	return x_out

# wrap atoms within simulation box (no jump option)
# input: x is coordinate sets of atoms
#		[ [x1, y1, z1], [x2, y2, z2], ...]
#	     box_xyz is box dimension
# output: a new wrapped coordinate inside box
def pbc_nojump_3d(xyz, box_xyz): 
	import numpy as np
	import copy
	xyz_out = copy.copy(xyz)
	#print(" io.pbc_nojump_1d:")
	for i_atom in range(len(xyz)):
		xyz_out[i_atom,0] = xyz[i_atom,0] - box_xyz[0]*np.floor(xyz[i_atom,0]/box_xyz[0]) # wrap coodinates within unit cell
		xyz_out[i_atom,1] = xyz[i_atom,1] - box_xyz[1]*np.floor(xyz[i_atom,1]/box_xyz[1]) # wrap coodinates within unit cell
		xyz_out[i_atom,2] = xyz[i_atom,2] - box_xyz[2]*np.floor(xyz[i_atom,2]/box_xyz[2]) # wrap coodinates within unit cell
	return xyz_out

# wrap atoms within simulation box (no jump option)
# test version because it is difficult to make all situations 
#  it would be slower than others because of comparison.
# input: x is coordinate sets of atoms
#		[ [x1, y1, z1], [x2, y2, z2], ...]
#	     box_xyz is box dimension
# output: a new wrapped coordinate inside box
def pbc_nojump_3d_test(xyz, box_xyz): 
	import numpy as np
	import copy
	xyz_out = copy.copy(xyz)
	temp = np.zeros(3)
	#print(" io.pbc_nojump_1d:")
	for i_atom in range(len(xyz)):
		for ix in range(3):
			temp[0] = xyz[i_atom,ix] - box_xyz[ix]*np.floor(xyz[i_atom,ix]/box_xyz[ix]) # wrap coodinates within unit cell
			temp[1] = xyz[i_atom,ix] - box_xyz[ix]*np.trunc(xyz[i_atom,ix]/box_xyz[ix]) # wrap coodinates within unit cell
			temp[2] = xyz[i_atom,ix] - box_xyz[ix]*np.ceil(xyz[i_atom,ix]/box_xyz[ix]) # wrap coodinates within unit cell
			xyz_out[i_atom,ix] = np.amin(np.absolute(temp))
	return xyz_out

# Boolean for Rectangular parallelepiped of unit cell
# input: unit_t is unit cell info (A of length, degree of angle) of atoms along time (t1, t2, ...)
#		[[x, y, z, angle_xy, angle_yz, angle_zx], [x, y, z, angle_xy, angle_yz, angle_zx], ...]
# output: true if Rectangular parallelepiped of unit cell. Otherwise, false.
# Example: check_rect_unit = rect_unit_cell(unit_t)

def rect_unit_cell(unit_t):
	check = True
	for unit in unit_t: 
		if check is False:
			break
		for i in range(3,5,1):	
			if check is False:
				break
			if unit[i] != 90.0:
				check = False
	return check

# Boolean for selected distance is beyond the half of the smallest box length
# input: unit_t is unit cell info (A of length, degree of angle) of atoms along time (t1, t2, ...)
#		[[x, y, z, angle_xy, angle_yz, angle_zx], [x, y, z, angle_xy, angle_yz, angle_zx], ...]
#			select_dist is what distance you select (A)
# output: true if select_dist is beyond the half. Otherwise, false.
# Example: check_beyond = dist_beyond_box(unit_t, select_dist)
# Assume: all angles are 90 deg.

def dist_beyond_box(unit_t, select_dist):
	check = False
	for unit in unit_t: 
		if check is True:
			break
		for i in range(0,2,1):	
			if check is True:
				break
			if select_dist >= unit[i]/2.0 :
				check = True
	return check

# get info of box_1d on axis
# input: box_1d is unit_cell trajectory on an axis 
#		[x_t1, x_t2, x_t3...]
# output: box_avg, box_std  
# Example: box_avg, box_std = box_1d(box_1d)
def box_1d(box_1d):
	print("coord.box_1d: ## delete soon ##")
	import numpy as np
	box_avg = np.mean(box_1d)
	box_std = np.std(box_1d)
	return box_avg, box_std

def box_1d_mode(box_1d, text, mode):
	print("coord.box_1d: ")
	import numpy as np
	box_avg = np.mean(box_1d)
	box_std = np.std(box_1d)
	if 'v' in mode:
		print(" {0} = {1:.5f} +- {2:.5f}".format(text,box_avg,box_std))
	return box_avg, box_std