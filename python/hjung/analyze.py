# make 1D histograms every single frame using 1d-coordinates
# input: x_t is 1d-coordinate of atoms along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
# 	 box_t is box dimension sets along time (t1, t2, ...)
#		[box_x(t0), box_x(t1), ...]
#        args.bin is which axis you want to do histogram
# output: histo_t is a new 1d-number profile trajectory
#	  bin_t is a bin position array for the 1d histogram 
# Example: number_t_1d, bin_t_1d = histo_t_1d(coordinates_1d, unit_cells_1d, args.bin) 

def histo_t_1d(x_t, box_x_t, bin_size):
	import numpy as np
	# check n_frames of x_t and box_x_t
	if len(x_t) != len(box_x_t):
		raise ValueError("# of time frame is not the same for input arrarys")
	# set number of bins (const) using input and initial box dimension
	n_bins = int(np.around(box_x_t[0]/bin_size))
	# initailize variables
	n_frames = len(x_t)
	histo_t = np.zeros((n_frames, n_bins))
	bin_t = np.zeros((n_frames, n_bins+1))
	# make histogram trajectory
	i_frame = 0
	for x in x_t:
		histo_t[i_frame], bin_t[i_frame] = np.histogram(x, bins=n_bins, range=(0,box_x_t[i_frame]))
		i_frame += 1
	return histo_t, bin_t

# make 1D histograms every single frame using 1d-coordinates with nbin
# input: x_t is 1d-coordinate of atoms along time (t1, t2, ...)
#			[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
# 		 box_t is box dimension sets along time (t1, t2, ...)
#			[box_x(t0), box_x(t1), ...]
#		 nbin is the number of bins you want 
# output: histo_t is a new 1d-number profile trajectory
#		  bin_t is a bin position array for the 1d histogram 
def histo_t_1d_nbin(x_t, box_x_t, nbin):
	import numpy as np
	# check n_frames of x_t and box_x_t
	if len(x_t) != len(box_x_t):
		raise ValueError("# of time frame is not the same for input arrarys")
	# set number of bins (const) using input and initial box dimension
	n_bins = nbin
	# initailize variables
	n_frames = len(x_t)
	histo_t = np.zeros((n_frames, n_bins))
	bin_t = np.zeros((n_frames, n_bins+1))
	# make histogram trajectory
	i_frame = 0
	for x in x_t:
		histo_t[i_frame], bin_t[i_frame] = np.histogram(x, bins=n_bins, range=(0,box_x_t[i_frame]))
		i_frame += 1
	return histo_t, bin_t

# make 1D histograms every single frame using 1d-coordinates with preset nbin
# input: x_t is 1d-coordinate of atoms along time (t1, t2, ...)
#			[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
# 		 box_t is box dimension sets along time (t1, t2, ...)
#			[box_x(t0), box_x(t1), ...]
#        bin_size is which axis you want to do histogram
#		 nbin is the number of bins you want 
#		 use_nbin is true if you want to use nbin instead of bin_size. otherwise, false
# output: histo_t is a new 1d-number profile trajectory
#		  bin_t is a bin position array for the 1d histogram 
def histo_t_1d_wbin_nbin(x_t, box_x_t, bin_size, nbin, use_nbin):
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
		histo_t[i_frame], bin_t[i_frame] = np.histogram(x, bins=n_bins, range=(0,box_x_t[i_frame]))
		i_frame += 1
	return histo_t, bin_t

# make 1d density histrogram of numbers
# input: x is 1d-array with numbers
# output: [[a collection set of density of each bins],[a collection set of bin]]]
# Example: prob_num_1d, prob_num_bin_1d = histo_num_dens_1d(number_1d)
def histo_num_dens_1d(x):
	import numpy as np
	n_bins = int(max(x))
	return np.histogram(x, bins=n_bins, density=True)

# get the information of interaction between two groups within cut-off distance at a given frame
# input: group1_coord or group2_coord is coordinates (position, A)for all atoms of group1 or group 2
#		[[x, y, z], [x, y, z], ...]
#			box is the unit cell info.
#		[x, y, z, angle_xy, angle_yz, angle_zx]
#			cutoff is cut-off distance (A)
# output: array of interaction within cut-off distance
#		[[group1_atom_x, group1_atom_y, group1_atom_z, unit_vector_x, unit_vector_y, unit_vector_z, 1/distance], ...] 
#		the sign of dist_xyz means the direction from group 1 to group 2, which is the same as sign of direction of coordinates
#		For example, if the an atom of group1 and group2 is left and right, the sign is positive (+). 
# Example: rdf_cut_array = rdf_cut(coordinates1[iframe], coordinates2[iframe], unit_cells2[iframe], args.cutoff)
# Assume: all information belongs in the same frame and the same box.
def rdf_cut_dt(group1_coord, group2_coord, box, cutoff):
	import numpy as np
	from numpy import linalg as LA

	# initial output array size (#possible combinations x 7 elements)
	output = np.zeros((len(group1_coord)*len(group2_coord), 7)) 
	#print('initializing rdf array')

	n_dist = 0 # total # selected distances
	i = 0
	for atom1 in group1_coord:
		if np.mod(i,50) == 0:
			print('...starting %d th atom in group1...' % i)
		for atom2 in group2_coord:
			dist = atom1 - atom2
			for xyz in range(0,2):
				dist[xyz] = dist[xyz] - np.around(dist[xyz]/box[xyz])*box[xyz]
			distance = LA.norm(dist)
			if distance <= cutoff:
				output[n_dist] = np.append(np.append(atom1,dist/distance),1.0/distance) # 7 elements
				n_dist += 1
		i += 1
	print('%d interactions are found' % n_dist)
	return output[0:n_dist]

# get the optimum block length for block average using tolerance of std of data
# all averages of blocks should be less than (average of total data +- tolerance*std)
# input: data_1d is 1D-data trajectory
#		[x_t1, x_t2, x_t3...]
#			tolerance is what you want to set acceptance line
# output: opt_block_length (integer) is the optimal block length. If not found, zero.  
# Example: opt_block_length = opt_block_length_1d(unit_cell_1d, tolerance)
# Assume: if mod(length of data, block_length) != 0, only last trajectory is considered, not beginning.
def opt_block_length_1d(data_1d, tolerance):
	print("Starting finding optimal block length....")
	import numpy as np
	ref_avg = np.mean(data_1d)
	ref_std = np.std(data_1d)
	ref_high = ref_avg + ref_std*tolerance
	ref_low = ref_avg - ref_std*tolerance
	print("Total data = %f +- %f \n" %(ref_avg,ref_std))
	# find max block length which have 5 blocks
	data_size = len(data_1d)
	quot, remain = divmod(data_size,5)
	max_block_length = quot
	opt_block_length = int(0)
	for block_length in range(2,max_block_length+1):
		quot, remain = divmod(data_size,block_length)
		last_data = data_1d[remain:]
		list_bavg = []
		for i in range(0,len(last_data),block_length):
			#  arithmetic mean of each block
			bavg = np.mean(last_data[i:i+block_length-1])
			list_bavg.append(bavg)
		# assume all blocks have the same weight
		pass_blength = True
		for j in range(0,len(list_bavg)):
			#print("list_bavg[%d] = %f, ref_high, %f, ref_low, %f" %(j, list_bavg[j], ref_high, ref_low))
			if list_bavg[j] > ref_high or list_bavg[j] < ref_low:
				#print("The value, %d, does not fit within %f +- %f.\n" 
				#	%(block_length, ref_avg, ref_std*tolerance))
				pass_blength = False
				break
			else:
				continue
		if pass_blength is False:
			continue
		else: # pass above check condition
			opt_block_length = block_length
			out_avg = np.mean(list_bavg)
			out_std = np.std(list_bavg)
			break
	if opt_block_length == 0:
		print("Not found optimal block length. \n")
		return 1
	else:
		print("Found optimal block length, %d." %opt_block_length)
		print("Before bloack average = %f +- %f " %(ref_avg,ref_std))
		print(" -> After bloack average = %f +- %f " %(out_avg,out_std))
		return opt_block_length


# get the optimum block length for block average using tolerance of std of data
# all averages of blocks should be less than (average of total data +- tolerance*std)
# input: data_t is array along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
#			block_length is block length you want
# output: new_data_t is averaged block array
# Example: number_t_1d = hjung.analyze.opt_block_length_1d(number_t_1d,block_length) 
# Assume: data_t matrix size is homogeneous. In other words, all element array in data_t is the same size
def block_average_1d(data_t, block_length):
	print("Starting block average....")
	import numpy as np
	nitems = len(data_t[0])
	nframes, remain = divmod(len(data_t),block_length)
	new_data_t = np.zeros((nframes,nitems)) # new data to store average
	temp = np.zeros(nitems) # temporary data set
	i = 0
	j = 0
	for iframe in range(remain,len(data_t)):
		temp = np.add(temp,data_t[iframe])
		i += 1
		if np.mod(i,block_length) == 0:
			new_data_t[j] = temp/float(block_length)
			j += 1 
			temp = np.zeros(nitems) # reset temp
	return new_data_t

