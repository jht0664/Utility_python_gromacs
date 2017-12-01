# merge all arraylists to one arraylist in order to remove time info.
# input: x_t is array along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
# output: [x1(t0), x2(t0), x3(t0), ..., x1(t1), x2(t1), x3(t1), ...] 
# Example: number_1d = merge_t_to_1d(number_t_1d)
def merge_t_to_1d(x_t):
	# internally merge arraylist in a list
	import itertools
	# merge 
	return [i for i in itertools.chain.from_iterable(x_t)]

# reduce unit cell array from MDAnalysis to 1d array along a given axis
#  Depending file extension, there are different order of unit cell component
# input: unitcells is array along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
# output: [x?(t0), x?(t1), x?(t2), x?(t3), ...] 
# Example: unit_cells_1d = reduce_unitcells_3d_to_1d(unit_cells, axis, structure, trajectory)
def reduce_unitcells_3d_to_1d(unit_cells, axis, structure, trajectory):
	import numpy as np
	print("="*30)
	print("reduce 3d unit cells to 1d: ######## delete soon #############")
	# gromace version
	if 'tpr' in structure and 'trr' in trajectory:
		print("We assume the input files are from Gromacs.")
		return unit_cells[:,axis]
	elif 'gro' in structure and 'trr' in trajectory:
		print("We assume the input files are from Gromacs.")
		return unit_cells[:,axis]
	# openmm version
	elif 'pdb' in structure and 'dcd' in input:
		print("We assume the input files are from OpenMM.")
		if axis == 0:
			return unit_cells[:,0]
		elif axis == 1:
			return unit_cells[:,2]
		elif axis == 2:
			return unit_cells[:,5]
		else:
			raise ValueError("wrong input of axis for histogram")
	# Monte Carlo version
	elif 'gro' in structure and 'xtc' in trajectory:
		print("We assume the input files are from Monte Carlo.")
		return unit_cells[:,axis]
	else:
		print("WARNING: unit_cells length may not be assigned correctly!")
		print("We use the same way for tpr, trr format)! You take risk.")
		return unit_cells[:,axis]	

def reduce_3d_to_1d(data_3d, axis):
	import numpy as np
	n_frames = len(data_3d)
	n_atoms = len(data_3d[0])
	data_1d = np.zeros((n_frames,n_atoms))
	n_frames = len(data_3d)
	for i_frame in range(len(data_3d)):
		for i_atom in range(n_atoms):
			#print("{} {} {}".format(i_frame,i_atom,axis))
			data_1d[i_frame][i_atom] = data_3d[i_frame][i_atom][axis]
	return data_1d

def convert_unitcell_3d(unit_cells, structure, trajectory):
	import numpy as np
	print("array.convert_unitcell_3d:")
	# gromace version
	if 'tpr' in structure and 'trr' in trajectory:
		print(" We assume the input files are from Gromacs.")
		return unit_cells
	elif 'gro' in structure and 'trr' in trajectory:
		print(" We assume the input files are from Gromacs.")
		return unit_cells
	# openmm version
	elif 'pdb' in structure and 'dcd' in input:
		print(" We assume the input files are from OpenMM.")
		output = np.zeros((len(unit_cells),3))
		output[:][0] = unit_cells[:,0]
		output[:][1] = unit_cells[:,2]
		output[:][2] = unit_cells[:,5]
		return output
	# Monte Carlo version
	elif 'gro' in structure and 'xtc' in trajectory:
		print(" We assume the input files are from Monte Carlo. If not, check if code works correctly")
		return unit_cells
	else:
		print(" WARNING: unit_cells length may not be assigned correctly!")
		print(" We use the same way for tpr, trr format)! You take risk.")
		return unit_cells