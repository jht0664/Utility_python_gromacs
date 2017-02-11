# Read 3D-coordinates from .trr trajectory files 
# input: tpr_filename, the filename of .tpr file in Gromacs
#	 trr_filename, the filename of .trr file in Gromacs
#	 select_atoms_filename, the filename including a command-line
#			for keeping trajectory of the selected atoms
# output: coordinates, xyz position of atoms, (#frame x #atoms x 3)
#	  unit_cells, box dimension (#frame x (x,y,z)) 
# Example: coordinates, unit_cells = read_coord_trr_3d('topol.tpr','traj.trr','b.select')

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
	if n_atoms == 0:
		raise ValueError("No atom is selected. %s may be wrong in grammer. Check again!" \
			% select_atoms_filename)
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

