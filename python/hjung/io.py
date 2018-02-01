# Read 3D-coordinates from .trr trajectory files 
# input: tpr_filename, the filename of .tpr file in Gromacs
#	 trr_filename, the filename of .trr file in Gromacs
#	 select_atoms_filename, the filename including a command-line
#			for keeping trajectory of the selected atoms
# output: coordinates, xyz position of atoms, (#frame x #atoms x 3)
#	  unit_cells, box dimension (#frame x (x,y,z)) 
# Example: coordinates, unit_cells = read_coord_trr_3d('topol.tpr','traj.trr','b.select')
def read_trr_3d_select1(tpr_filename, trr_filename, select_atoms_filename, mode):
	print("read_trr_3d_select1:")
	# import
	import MDAnalysis
	import numpy as np
	# check the arg, mode
	outmode = np.zeros(3,dtype=bool) # use true and false in python as 1 and 0
	if 'pos' in mode:
		outmode[0] = True
	if 'vel' in mode:
		outmode[1] = True
	if 'forc' in mode:
		outmode[2] = True
	if not np.any(outmode): # if all is false
		raise ValueError(" wrong arg mode {}".format(mode))
	ndata = sum(bool(x) for x in outmode)
	print(" output data #sets = {} by your mode setting {} ".format(ndata,mode))
	# read a line of select command-line for MDAnalysis
	select_command = []
	if select_atoms_filename is not None:
		try:
			open_file = open(select_atoms_filename, 'r')
		except IOError:
			raise IOError(" problem with opening ",select_atoms_filename)
		select_command_temp = open_file.readline().strip()
		open_file.close()
		print(" select written in {}: {}".format(select_atoms_filename,select_command_temp))
		select_command.append(select_command_temp)
	else:
		raise ValueError(" wrong select atom files {}".format(select_atoms_filename))
	# Read trajectory using MDAnalysis 
	u = MDAnalysis.Universe(tpr_filename,trr_filename)
	n_frames = len(u.trajectory)
	# obtain a set of atom index
	n_atoms = []
	atoms = []
	for iselect in select_command: 	
		list_atoms = u.select_atoms(iselect).indices
		if len(list_atoms) == 0:
			raise ValueError(" No atom is selected. {} may be wrong in grammer.".format(iselect))
		atoms.append(list_atoms)
		n_atoms.append(len(list_atoms))
	print(" selected total #atoms: {}".format(n_atoms))	
	# initailize variables
	if ndata > 1:
		data1 = np.zeros((ndata, n_frames, n_atoms[0], 3))
	else:
		data1 = np.zeros((n_frames, n_atoms[0], 3))
	unit_cells = np.zeros((n_frames, 6))
	# read trajectory
	print(" starting reading trajectory...")
	i_frame = 0
	mod_frame = process_init()
	for ts in u.trajectory:
		try:
			if ndata == 1:
				if outmode[0]:
					tmp = np.array(ts._pos)
				if outmode[1]:
					tmp = np.array(ts._velocities)
				if outmode[2]:	
					tmp = np.array(ts._forces)
				data1[i_frame, :, :] = tmp[atoms[0]]
			else:
				dataset = 0
				if outmode[0]:
					tmp = np.array(ts._pos)
					data1[dataset, i_frame, :, :] = tmp[atoms[0]]
					dataset = dataset + 1
				if outmode[1]:
					tmp = np.array(ts._velocities)
					data1[dataset, i_frame, :, :] = tmp[atoms[0]]
					dataset = dataset + 1
				if outmode[2]:	
					tmp = np.array(ts._forces)
					data1[dataset, i_frame, :, :] = tmp[atoms[0]]
					dataset = dataset + 1
				if dataset != ndata:
					raise ValueError(" weird number of data set reading trajectory. {} {}".format(dataset,ndata))
			unit_cells[i_frame, :] = ts._unitcell
		except IndexError:
			raise ValueError(" There are more coordinates to be read than indicated in the header.")
		i_frame += 1
		mod_frame = process_print(i_frame,n_frames,mod_frame)
	# check consistency; final i_frame should be the same as # frames
	if i_frame != n_frames:
		print(" actual nframes {} in trajectory != the length claimed in header of trajectory {}".format(i_frame, n_frames))
		print(" saving trajectory is problem (due to limit of disk quota). Size of your data will be {} by force".format(i_frame))
	print("# frames will be extracted from trajectory = {} (excludes t=0)".format(i_frame-1))
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
	if ndata == 1:
		return data1[1:i_frame], unit_cells[1:i_frame] #return data1[0:i_frame-1], data2[0:i_frame-1], unit_cells[0:i_frame-1]
	else:
		return data1[:,1:i_frame,:,:], unit_cells[1:i_frame] # return data1[0]:pos, data1[1]:force if your mode is "pos force"

# Read 2 3D-coordinates from .trr trajectory files 
# input: tpr_filename, the filename of .tpr file in Gromacs
#	 	trr_filename, the filename of .trr file in Gromacs
#	 	select_atoms_filename#, the filename including a command-line
#			for keeping trajectory of the selected atoms
# output: coordinates#, xyz position of atoms, (#frame x #atoms x 3)
#	  	unit_cells, box dimension (#frame x (x,y,z)) 
# Example: coordinates1, coordinates2, unit_cells = read_coord_trr_3d_select2('topol.tpr','traj.trr','b.select','peo.select')
def read_trr_3d_select2(tpr_filename, trr_filename, select_atoms_filename1, select_atoms_filename2, mode):
	print("io.read_trr_3d_select2:")
	# import
	import MDAnalysis
	import numpy as np
	# check the arg, mode
	outmode = np.zeros(3,dtype=bool) # use true and false in python as 1 and 0
	if 'pos' in mode:
		outmode[0] = True
	if 'vel' in mode:
		outmode[1] = True
	if 'forc' in mode:
		outmode[2] = True
	if not np.any(outmode): # if all is false
		raise ValueError(" wrong arg mode {}".format(mode))
	ndata = sum(bool(x) for x in outmode)
	print(" output data #sets = {} by your mode setting {} ".format(ndata,mode))
	# read a line of select command-line for MDAnalysis
	select_command = []
	for select_atoms_filename in [select_atoms_filename1, select_atoms_filename2]:
		if select_atoms_filename is not None:
			try:
				open_file = open(select_atoms_filename, 'r')
			except IOError:
				raise IOError(" problem with opening ",select_atoms_filename)
			select_command_temp = open_file.readline().strip()
			open_file.close()
			print(" select written in {}: {}".format(select_atoms_filename,select_command_temp))
			select_command.append(select_command_temp)
		else:
			raise ValueError(" wrong select atom files {}".format(select_atoms_filename))
	# Read trajectory using MDAnalysis 
	u = MDAnalysis.Universe(tpr_filename,trr_filename)
	n_frames = len(u.trajectory)
	# obtain a set of atom index
	n_atoms = []
	atoms = []
	for iselect in select_command: 	
		list_atoms = u.select_atoms(iselect).indices
		print(" You selected {} atoms for selection {}".format(len(list_atoms),iselect))
		if len(list_atoms) == 0:
			raise ValueError(" No atom is selected. {} may be wrong in grammer.".format(iselect))
		atoms.append(list_atoms)
		n_atoms.append(len(list_atoms))
	print(" selected total #atoms: {}".format(n_atoms))	
	# initailize variables
	if ndata > 1:
		data1 = np.zeros((ndata, n_frames, n_atoms[0], 3))
		data2 = np.zeros((ndata, n_frames, n_atoms[1], 3))
	else:
		data1 = np.zeros((n_frames, n_atoms[0], 3))
		data2 = np.zeros((n_frames, n_atoms[1], 3))
	unit_cells = np.zeros((n_frames, 6))
	# read trajectory
	print(" starting reading trajectory...")
	i_frame = 0
	mod_frame = process_init()
	for ts in u.trajectory:
		try:
			if ndata == 1:
				if outmode[0]:
					tmp = np.array(ts._pos)
				if outmode[1]:
					tmp = np.array(ts._velocities)
				if outmode[2]:	
					tmp = np.array(ts._forces)
				data1[i_frame, :, :] = tmp[atoms[0]]
				data2[i_frame, :, :] = tmp[atoms[1]]
			else:
				dataset = 0
				if outmode[0]:
					tmp = np.array(ts._pos)
					data1[dataset, i_frame, :, :] = tmp[atoms[0]]
					data2[dataset, i_frame, :, :] = tmp[atoms[1]]
					dataset = dataset + 1
				if outmode[1]:
					tmp = np.array(ts._velocities)
					data1[dataset, i_frame, :, :] = tmp[atoms[0]]
					data2[dataset, i_frame, :, :] = tmp[atoms[1]]
					dataset = dataset + 1
				if outmode[2]:	
					tmp = np.array(ts._forces)
					data1[dataset, i_frame, :, :] = tmp[atoms[0]]
					data2[dataset, i_frame, :, :] = tmp[atoms[1]]
					dataset = dataset + 1
				if dataset != ndata:
					raise ValueError(" weird number of data set reading trajectory. {} {}".format(dataset,ndata))
			unit_cells[i_frame, :] = ts._unitcell
		except IndexError:
			raise ValueError(" There are more coordinates to be read than indicated in the header.")
		i_frame += 1
		mod_frame = process_print(i_frame,n_frames,mod_frame)
	# check consistency; final i_frame should be the same as # frames
	if i_frame != n_frames:
		print(" actual nframes {} in trajectory != the length claimed in header of trajectory {}".format(i_frame, n_frames))
		print(" saving trajectory is problem (due to limit of disk quota). Size of your data will be {} by force".format(i_frame))
	print(" # frames will be extracted from trajectory = {} (excludes t=0)".format(i_frame-1))
    # box info 
	if all(unit_cells[0,:] == unit_cells[1,:]):
		print(" The system may be in NVT ensemble")
	else:
		# for gromacs (tpr, trr files)
		# unit_cells = [length_x, length_y, length_z, angles, ...]
		if 'trr' in trr_filename and 'tpr' in tpr_filename:
			if unit_cells[0][0] == unit_cells[1][0] and unit_cells[0][1] == unit_cells[1][1]:
				print(" may be in NPAT ensemble")
			else:
				print(" may be in NPT ensemble")
		# for openmm (pdb, dcd files)
		# unit_cells = [length_x, alpha angle, length_y, beta angle, theta angle, length_z]
		if 'dcd' in trr_filename and 'pdb' in tpr_filename:
			if unit_cells[0][0] == unit_cells[1][0] and unit_cells[0][2] == unit_cells[1][2]:
				print(" may be in NPAT ensemble")
			else:
				print(" may be in NPT ensemble")
	if ndata == 1:
		#return data1[0:i_frame-1], data2[0:i_frame-1], unit_cells[0:i_frame-1]
		return data1[1:i_frame], data2[1:i_frame], unit_cells[1:i_frame]
	else:
		return data1[:,1:i_frame,:,:], data2[:,1:i_frame,:,:], unit_cells[1:i_frame]		

# Read some 3D-coordinates from .trr trajectory files 
# input: tpr_filename, the filename of .tpr file in Gromacs
#	 	 trr_filename, the filename of .trr file in Gromacs
#	 	 select_atoms_filename#, the filename including a command-line
#			for keeping trajectory of the selected atoms
# output: coordinates#, xyz position of atoms, (#frame x #atoms x 3)
#	  	  unit_cells, box dimension (#frame x (x,y,z)) 
# Example: coordinates, unit_cells = read_coord_trr_3d_select2('topol.tpr','traj.trr','b.select','peo.select')
def read_trr_3d_selects(tpr_filename, trr_filename, select_filename, mode):
	print("io.read_trr_3d_selects: ######## INCOMPLETE coding #### ")
	# import
	import MDAnalysis
	import numpy as np
	# check the arg, mode
	outmode = np.zeros(3,dtype=bool) # use true and false in python as 1 and 0
	if 'pos' in mode:
		outmode[0] = True
	if 'vel' in mode:
		outmode[1] = True
	if 'forc' in mode:
		outmode[2] = True
	if not np.any(outmode): # if all is false
		raise ValueError(" wrong arg mode {}".format(mode))
	ndata = sum(bool(x) for x in outmode)
	print(" output data #sets = {} by your mode setting {} ".format(ndata,mode))
	# read a line of select command-line for MDAnalysis
	select_command = []
	if select_filename is not None:
		try:
			open_file = open(select_filename, 'r')
		except IOError:
			raise IOError(" problem with opening ",select_atoms_filename)
		for line in oepn_file:
			iline = line.strip()
			# line is blank or comment line
			if not iline or iline.startswith('#') or iline.startswith('@'):
				continue
			print(" select written in {}: {}".format(select_filename,iline))
			select_command.append(iline)
		n_select = len(select_command)
		open_file.close()
	else:
		raise ValueError(" wrong select atom files {}".format(select_atoms_filename))
	# Read trajectory using MDAnalysis 
	u = MDAnalysis.Universe(tpr_filename,trr_filename)
	n_frames = len(u.trajectory)
	# obtain a set of atom index
	n_atoms = []
	atoms = []
	for iselect in select_command: 	
		list_atoms = u.select_atoms(iselect).indices
		if len(list_atoms) == 0:
			raise ValueError(" No atom is selected. {} may be wrong in grammer.".format(iselect))
		atoms.append(list_atoms)
		n_atoms.append(len(list_atoms))
	print(" selected total #atoms: {}".format(n_atoms))	
	######################################## followings are not complete
	# initailize variables with zeros
	if ndata > 1:
		data1 = np.zeros((ndata, n_frames, n_atoms[0], 3))
		data2 = np.zeros((ndata, n_frames, n_atoms[1], 3))
	else:
		data1 = np.zeros((n_frames, n_atoms[0], 3))
		data2 = np.zeros((n_frames, n_atoms[1], 3))
	unit_cells = np.zeros((n_frames, 6))
	# read trajectory
	print(" starting reading trajectory...")
	i_frame = 0
	mod_frame = process_init()
	for ts in u.trajectory:
		try:
			if ndata == 1:
				if outmode[0]:
					tmp = np.array(ts._pos)
				if outmode[1]:
					tmp = np.array(ts._velocities)
				if outmode[2]:	
					tmp = np.array(ts._forces)
				data1[i_frame, :, :] = tmp[atoms[0]]
				data2[i_frame, :, :] = tmp[atoms[1]]
			else:
				dataset = 0
				if outmode[0]:
					tmp = np.array(ts._pos)
					data1[dataset, i_frame, :, :] = tmp[atoms[0]]
					data2[dataset, i_frame, :, :] = tmp[atoms[1]]
					dataset = dataset + 1
				if outmode[1]:
					tmp = np.array(ts._velocities)
					data1[dataset, i_frame, :, :] = tmp[atoms[0]]
					data2[dataset, i_frame, :, :] = tmp[atoms[1]]
					dataset = dataset + 1
				if outmode[2]:	
					tmp = np.array(ts._forces)
					data1[dataset, i_frame, :, :] = tmp[atoms[0]]
					data2[dataset, i_frame, :, :] = tmp[atoms[1]]
					dataset = dataset + 1
				if dataset != ndata:
					raise ValueError(" weird number of data set reading trajectory. {} {}".format(dataset,ndata))
			unit_cells[i_frame, :] = ts._unitcell
		except IndexError:
			raise ValueError(" There are more coordinates to be read than indicated in the header.")
		i_frame += 1
		mod_frame = process_print(i_frame,n_frames,mod_frame)
	# check consistency; final i_frame should be the same as # frames
	if i_frame != n_frames:
		print(" actual nframes {} in trajectory != the length claimed in header of trajectory {}".format(i_frame, n_frames))
		print(" saving trajectory is problem (due to limit of disk quota). Size of your data will be {} by force".format(i_frame))
	print(" # frames will be extracted from trajectory = {} (excludes t=0)".format(i_frame-1))
    # box info 
	if all(unit_cells[0,:] == unit_cells[1,:]):
		print(" The system may be in NVT ensemble")
	else:
		# for gromacs (tpr, trr files)
		# unit_cells = [length_x, length_y, length_z, angles, ...]
		if 'trr' in trr_filename and 'tpr' in tpr_filename:
			if unit_cells[0][0] == unit_cells[1][0] and unit_cells[0][1] == unit_cells[1][1]:
				print(" may be in NPAT ensemble")
			else:
				print(" may be in NPT ensemble")
		# for openmm (pdb, dcd files)
		# unit_cells = [length_x, alpha angle, length_y, beta angle, theta angle, length_z]
		if 'dcd' in trr_filename and 'pdb' in tpr_filename:
			if unit_cells[0][0] == unit_cells[1][0] and unit_cells[0][2] == unit_cells[1][2]:
				print(" may be in NPAT ensemble")
			else:
				print(" may be in NPT ensemble")
	if ndata == 1:
		#return data1[0:i_frame-1], data2[0:i_frame-1], unit_cells[0:i_frame-1]
		return data1[1:i_frame], data2[1:i_frame], unit_cells[1:i_frame]
	else:
		return data1[:,1:i_frame,:,:], data2[:,1:i_frame,:,:], unit_cells[1:i_frame]		

# rename the existing filename if the filename already exists
# input: filename is a filename what you want to rename if the filename already exists 
# Example: rename_existing_file("../../fit.plot")
# ref: http://code.activestate.com/recipes/578116-move-files-with-rename-if-required/
def rename_existing_file(filename):
	import sys, os
	print('backup_existing_file:')

	if not os.path.exists(filename):
		print(" %s does not exist." %filename)
		return
	# get info
	dirname = os.path.dirname(filename)
	basename = os.path.basename(filename) # basename is a fiilname excluding path. e.g. ../../fit.plot -> fit.plot
	head, tail = os.path.splitext(basename) # head = fit, tail = .plot
	#dst_file = os.path.join(dst_dir, basename) # dst_dir is target folder path

	# find better filename (= filename2)
	filename2 = filename
	count = 0
	while os.path.exists(filename2):
		filename2 = os.path.join(dirname, '{0}.{1}{2}'.format(head, count, tail))
		count += 1
	print(" rename {0} to {1}".format(filename,filename2))
	os.rename(filename, filename2)


# initialize processing bar
# output: process tiem and wall time
# Example: mod_frame = process_init()
def process_init():
	print("io.process_init: ")
	return 10 # initial mode number = 10

# print process bar
# input: itime, current time
#        ftime, final time
#        mod_time, frequency of printing.
# output: mod_time, new or old printing frequency 
# Example: mod_frame = process_print(itime+1, ftime, mod_frame)
def process_print(itime, ftime, mod_time):
	if itime%mod_time == 0:
		print("... {0} th frame reading ({1:.0%}) ...".format(itime,itime/ftime))
		if (itime/mod_time)%10 == 0:
			mod_time = mod_time*10
	return mod_time

# read mass values
# input: filename
# output: mw, molecular weight
#         divider, number of selected atoms per molecule or polymer chain length
# Example: mw, divider = read_mass2(filename)
def read_mass2(filename):
	print("io.read_mass2:")
	import numpy as np
	# read mass file for two selection
	print(" reading the mass file, {}".format(filename))
	try:
		massinfo = open(filename, 'r')
	except IOError:
		print(" Problem with opening ",filename)
		exit()
	divider = []
	mw = []
	for line in massinfo:
		line = line.strip()
		line_m = line.rsplit()
		divider.append(float(line_m[0]))
		mw.append(float(line_m[1]))
	massinfo.close()
	# type
	divider = np.array(divider,dtype=np.float)
	mw = np.array(mw,dtype=np.float)
	if len(divider) != 2 or len(mw) != 2:
		ValueError("Wrong format in {} file".format(filename))
	print(" dividers[select1,select2] = {}".format(divider))
	print(" mw[select1,select2] = {}".format(mw))
	return mw, divider