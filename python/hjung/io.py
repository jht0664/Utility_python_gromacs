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
	nmod = 10
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
		if i_frame%nmod == 0:
			print("... {0} th frame reading ({1:.0%}) ...".format(i_frame,i_frame/n_frames))
			if (i_frame/nmod)%10 == 0:
				nmod = nmod*10
	# check consistency; final i_frame should be the same as # frames
	if i_frame != n_frames:
		raise ValueError("# of frames to read %d does not agree with the length of trajectory file %d" \
	                             % (i_frame, n_frames))
	# box info 
	if all(unit_cells[0,:] == unit_cells[1,:]):
		print("The system may be in NVT ensemble")
	else:
		# for gromacs (tpr, trr files)
		# unit_cells = [length_x, length_y, length_z, angles, ...]
		if 'trr' in trr_filename and 'tpr' in tpr_filename:
			if unit_cells[0][0] == unit_cells[1][0] and unit_cells[0][1] == unit_cells[1][1]:
				print("may be in NPAT ensemble")
			else:
				print("may be in NPT ensemble")
		# for openmm (pdb, dcd files)
		# unit_cells = [length_x, alpha angle, length_y, beta angle, theta angle, length_z]
		if 'dcd' in trr_filename and 'pdb' in tpr_filename:
			if unit_cells[0][0] == unit_cells[1][0] and unit_cells[0][2] == unit_cells[1][2]:
				print("may be in NPAT ensemble")
			else:
				print("may be in NPT ensemble")
	return coordinates, unit_cells


# Read 2 3D-coordinates from .trr trajectory files 
# input: tpr_filename, the filename of .tpr file in Gromacs
#	 	trr_filename, the filename of .trr file in Gromacs
#	 	select_atoms_filename#, the filename including a command-line
#			for keeping trajectory of the selected atoms
# output: coordinates#, xyz position of atoms, (#frame x #atoms x 3)
#	  	unit_cells, box dimension (#frame x (x,y,z)) 
# Example: coordinates1, coordinates2, unit_cells = read_coord_trr_3d_select2('topol.tpr','traj.trr','b.select','peo.select')

def read_coord_trr_3d_select2(tpr_filename, trr_filename, select_atoms_filename1, select_atoms_filename2):
	# print function name
	print("read_coord_trr_3d_select2:")
	# read a line of select command-line in MDAnalysis
	select_command = []
	for select_atoms_filename in [select_atoms_filename1, select_atoms_filename2]:
		if select_atoms_filename is not None:
			try:
				open_file = open(select_atoms_filename, 'r')
			except IOError:
				raise IOError("Problem with opening ",select_atoms_filename)
			print("reading the first line in %s file" %select_atoms_filename)
			select_command_temp = open_file.readline().strip()
			print("select: %s" %select_command_temp)
			select_command.append(select_command_temp)
			open_file.close()

	# import
	import MDAnalysis
	import numpy as np
	
	# Read trajectory using MDAnalysis 
	u = MDAnalysis.Universe(tpr_filename,trr_filename)
	n_frames = len(u.trajectory)
	print("# frames = ", n_frames)

	# selected atoms
	n_atoms = []
	atoms = []
	for iselect in select_command: 	
		list_atoms = u.select_atoms(iselect).indices
		if len(list_atoms) == 0:
			raise ValueError("No atom is selected. %s may be wrong in grammer. Check again!" \
				% iselect)
		atoms.append(list_atoms)
		n_atoms.append(len(list_atoms))
	print("# atoms: %s" %n_atoms)	
	
	# initailize variables
	coordinates1 = np.zeros((n_frames, n_atoms[0], 3))
	coordinates2 = np.zeros((n_frames, n_atoms[1], 3))
	#velocities = np.zeros(n_frames, n_atoms, 3)
	#forces = np.zeros(n_frames, n_atoms, 3)
	unit_cells = np.zeros((n_frames, 6))
	
	# read trajectory
	print("Starting reading trajectory...")
	i_frame = 0
	nmod = 10
	for ts in u.trajectory:
		try:
			tmp = np.array(ts._pos)
			coordinates1[i_frame, :, :] = tmp[atoms[0]] # read positions in 10^-10 m (A)
			coordinates2[i_frame, :, :] = tmp[atoms[1]] # read positions in 10^-10 m (A)
			#np.append(velocities, ts._velocities)
			#np.append(forces, ts._forces)
			unit_cells[i_frame, :] = ts._unitcell
		except IndexError:
			raise ValueError("There are more coordinates to be read than indicated in the header.")
		i_frame += 1
		if i_frame%nmod == 0:
			print("... {0} th frame reading ({1:.0%}) ...".format(i_frame,i_frame/n_frames))
			if (i_frame/nmod)%10 == 0:
				nmod = nmod*10
	# check consistency; final i_frame should be the same as # frames
	if i_frame != n_frames:
		raise ValueError("# of frames to read %d does not agree with the length of trajectory file %d" \
	                             % (i_frame, n_frames))

	# box info 
	if all(unit_cells[0,:] == unit_cells[1,:]):
		print("The system may be in NVT ensemble")
	else:
		# for gromacs (tpr, trr files)
		# unit_cells = [length_x, length_y, length_z, angles, ...]
		if 'trr' in trr_filename and 'tpr' in tpr_filename:
			if unit_cells[0][0] == unit_cells[1][0] and unit_cells[0][1] == unit_cells[1][1]:
				print("may be in NPAT ensemble")
			else:
				print("may be in NPT ensemble")
		# for openmm (pdb, dcd files)
		# unit_cells = [length_x, alpha angle, length_y, beta angle, theta angle, length_z]
		if 'dcd' in trr_filename and 'pdb' in tpr_filename:
			if unit_cells[0][0] == unit_cells[1][0] and unit_cells[0][2] == unit_cells[1][2]:
				print("may be in NPAT ensemble")
			else:
				print("may be in NPT ensemble")
				
	return coordinates1, coordinates2, unit_cells
