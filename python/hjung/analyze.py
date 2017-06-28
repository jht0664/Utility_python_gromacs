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

# The best way to get all the divisors of a number
# http://stackoverflow.com/questions/171765/what-is-the-best-way-to-get-all-the-divisors-of-a-number
def divisors(number):
	n = 1
	list_n = []
	while(n<number):
		if(number%n==0):
			list_n.append(n)
		else:
			pass
		n += 1
	return list_n

# get the optimum block length for block average using tolerance of std of data
# all averages of blocks should be less than (average of total data +- tolerance*std)
# input: data_1d is 1D-data trajectory
#		[x_t1, x_t2, x_t3...]
#			tolerance is what you want to set acceptance line
# output: opt_block_length (integer) is the optimal block length. If not found, zero.  
# Example: opt_block_length = opt_block_length_1d(unit_cell_1d, tolerance)
# Assume: if mod(length of data, block_length) != 0, only last trajectory is considered, not beginning.
def opt_block_length_1d(data_1d, tolerance):
	import numpy as np
	print("# will be removed later")

	#print("opt_block_length_1d: ")
	
	# check tolerance	
	if tolerance > 1: 
		return int(tolerance)
	elif tolerance < 0:
		return 1

	# info before block average
	ref_avg = np.mean(data_1d)
	ref_std = np.std(data_1d)
	ref_high = ref_avg + ref_std*tolerance
	ref_low = ref_avg - ref_std*tolerance
	#print("(before block avg) Total data = %f +- %f" %(ref_avg,ref_std))
	if (100*float(ref_std)/ref_avg) < 0.05: # std is very small (<0.05%)
	#	print("original std value is less than 5%. block average not necessary.\n")
		return 2
	
	# make list of possible block lengths
	data_size = len(data_1d)
	quot, remain = divmod(data_size,2)
	if remain != 0: # if odd frame, force to make even frame
		data_size -= 1
	poss_block_length = divisors(data_size) # make all list of divisors of #frames
	poss_block_length = poss_block_length[1:-3] # remove 1 and itself
	if len(poss_block_length) == 0: # if the list is small, return no average
		return 2

	# determine optimal block length using the list
	opt_block_length = int(0)
	for block_length in poss_block_length:
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
	#	print("Not found optimal block length.\n")
		return 2
	else:
		#print("Found optimal block length, %d." %opt_block_length)
		#print("Before bloack average = %f +- %f " %(ref_avg,ref_std))
		#print(" -> After bloack average = %f +- %f " %(out_avg,out_std))
		return opt_block_length

# get the optimum block length for block average using tolerance of std of data
# all averages of blocks should be less than (average of total data +- tolerance*std)
# input: data_1d is 1D-data trajectory
#		[x_t1, x_t2, x_t3...]
#			tolerance is what you want to set acceptance line
# output: opt_block_length (integer) is the optimal block length. If not found, zero.  
# Example: opt_block_length = opt_block_length_1d(unit_cell_1d, tolerance)
# Assume: if mod(length of data, block_length) != 0, only last trajectory is considered, not beginning.
def opt_block_length_1d_t(data_1d_t, tolerance):
	print("opt_block_length_1d_t: # will be removed later")
	import numpy as np
	if tolerance > 0.0 and tolerance <= 1.0:
		print("optimize block length under tolerance")
		data_t_1d = data_1d_t.transpose()
		opt_lengths = []
		for data_t in data_t_1d:
			temp = opt_block_length_1d(data_t,tolerance) 
			opt_lengths.append(temp)
		print("Here we found this length is optimum %d" %(np.amax(opt_lengths)))
		return np.amax(opt_lengths)
	elif tolerance > 1: 
		return int(tolerance)
	elif tolerance == 0.0:
		return 1 # No Block Average
	else:
		raise ValueError("Wrong tolerance, %s, for block length" %tolerance)

# get the optimum block length for block average using tolerance of std of data
# all averages of blocks should be less than (average of total data +- tolerance*std)
# (x1(t0)+x1(t1)+x1(t2)+....)/block_length = block average
# input: data_t is array along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
#			block_length is block length you want
# output: new_data_t is averaged block array
# Example: number_t_1d = hjung.analyze.opt_block_length_1d(number_t_1d,block_length) 
# Assume: data_t matrix size is homogeneous. In other words, all element array in data_t is the same size
def block_average_1d(data_t, block_length):
	import numpy as np
	print("block_average_1d: # will be removed later")
	
	if block_length == 1: # no block average
		print(" No block average")
		return data_t
	
	nitems = len(data_t[0])
	nframes, remain = divmod(len(data_t),block_length)
	print("set block length, %d" %block_length)
	print("uses %d blocks" %nframes)
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

# 1D auto-correlation function 
# \Delta (\Delta N(i,t))\quad =\quad \Delta N(i,t)\quad -\quad { <\Delta N(i,t)> }_{ i }
# input: data_1d_t is array along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ..., xn(t0)], [x1(t1), x2(t1), x3(t1), ..., xn(t1)], ...]
#		 setmode is setting for ndimage.correlate(mode). See below table.
# output: [[tau_x1(t0), tau_x2(t0), tau_x3(t0), ..., tau_xn(t0)], [tau_x1(t1), tau_x2(t1), tau_x3(t1),..., tau_xn(t1)], ...]
#          same size. But different meaning of x's. tau_x range [-half_total_t:half_total_t]
# Example: autocorr_1d = hjung.analyze.autocorr_1d_t(data_1d_t,mode) 
# Assume: data_t matrix size is homogeneous. The value after xn(t0) is x1(t0). 
# See detail ndimage.correlate(mode='wrap')
################################################################
##  mode       |   Ext   |         Input          |   Ext
##  -----------+---------+------------------------+---------
##  'mirror'   | 4  3  2 | 1  2  3  4  5  6  7  8 | 7  6  5
##  'reflect'  | 3  2  1 | 1  2  3  4  5  6  7  8 | 8  7  6
##  'nearest'  | 1  1  1 | 1  2  3  4  5  6  7  8 | 8  8  8
##  'constant' | 0  0  0 | 1  2  3  4  5  6  7  8 | 0  0  0
##  'wrap'     | 6  7  8 | 1  2  3  4  5  6  7  8 | 1  2  3
################################################################
def autocorr_1d_t(data_1d_t,setmode):
	print("Starting Autocorrelation 1d with mode %s ...." %setmode)
	
	import numpy as np
	from scipy import ndimage

	acf_data = []
	for data_1d in data_1d_t:
		delta_data_1d = data_1d - data_1d.mean()
		acf_data_1d = ndimage.correlate(delta_data_1d,delta_data_1d,mode=setmode)
		acf_data_1d /= (data_1d.var()*len(delta_data_1d)) # normalize
		acf_data.append(acf_data_1d)	
	return acf_data

# 1D auto-correlation function using fft
def autocorr_1d_fft(data_1d):
	# Compute autocorrelation using FFT
	# The idea comes from 
	# http://dsp.stackexchange.com/a/1923/4363 (Hilmar)
	import numpy as np
	fft = np.fft
	x = np.asarray(x)
	N = len(x)
	x = x-x.mean()
	s = fft.fft(x, N*2-1)
	result = np.real(fft.ifft(s * np.conjugate(s), N*2-1))
	result = result[:N]
	result /= result[0]
	return result

# 1D auto-correlation function using fft with data_1d_t
def autocorr_1d_fft_t(data_1d_t):
	acf_data = []
	for data_1d in data_1d_t:
		acf_data_1d = autocorr_1d_fft(data_1d)
		acf_data.append(acf_data_1d)
	return acf_data

# 1D convolution with data#_t which are the same size
# input: data#_1d_t is array along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ..., xn(t0)], [x1(t1), x2(t1), x3(t1), ..., xn(t1)], ...]
#       data1 is ideal case or reference and data2 is what you are going to fit into data1.
#		 setmode is setting for ndimage.convolve(mode). See below table.
# output: [max_value_index(t0), max_value_index(t1), ...]
#          returns the index of max_value at each frame 
# Example: align_shift = hjung.analyze.convolve_1d_t_max(data1_1d_t,data2_1d_t,setmode) 
# See detail about modes
################################################################
##  mode       |   Ext   |         Input          |   Ext
##  -----------+---------+------------------------+---------
##  'mirror'   | 4  3  2 | 1  2  3  4  5  6  7  8 | 7  6  5
##  'reflect'  | 3  2  1 | 1  2  3  4  5  6  7  8 | 8  7  6
##  'nearest'  | 1  1  1 | 1  2  3  4  5  6  7  8 | 8  8  8
##  'constant' | 0  0  0 | 1  2  3  4  5  6  7  8 | 0  0  0
##  'wrap'     | 6  7  8 | 1  2  3  4  5  6  7  8 | 1  2  3
################################################################
def convolve_1d_t_min(data1_1d_t, data2_1d_t, setmode):
	import numpy as np
	from scipy import ndimage

	# check total # element and length of datas
	data1_1d_t = np.array(data1_1d_t)
	data2_1d_t = np.array(data2_1d_t)
	if data1_1d_t.size != data2_1d_t.size:
		raise ValueError("Error: # elements of datas are not same.")
	len_data1 = len(data1_1d_t)
	len_data2 = len(data2_1d_t)
	if len_data1 != len_data2:
		raise ValueError("Error: length of datas are not same.")
	if len_data1*len(data1_1d_t[0]) != data1_1d_t.size \
		or len_data2*len(data2_1d_t[0]) != data2_1d_t.size:
		raise ValueError("Error: datas are not homogeneous shape.")
	
	# Do convolution at the same frame
	convolve_data = []
	for iframe in range(len_data1):
		convolve_temp = ndimage.convolve(data1_1d_t[iframe],data2_1d_t[iframe],mode=setmode)
		# No normalization
		convolve_data.append(convolve_temp)

	return np.argmin(convolve_data, axis=1)

# Align data_1d_t using convolution with (spatial) autocorrelation function of data_1d_t
# input: data_1d_t is array along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ..., xn(t0)], [x1(t1), x2(t1), x3(t1), ..., xn(t1)], ...]
#		is what you are going to fit into autocorreflation function.
#		 acf_1d_t is array along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ..., xn(t0)], [x1(t1), x2(t1), x3(t1), ..., xn(t1)], ...]
#       is (spatial) autocorreflation function.
#		setmode is setting for ndimage.convolve(mode). See below table.
# output: [[x3(t0), x4(t0), x5(t0), ..., x1(t0)], [x8(t1), x9(t1), x10(t1), ..., x1(t1)], ...]
#          returns new aligned data_1d_t
# Example: new_data_1d_t = hjung.analyze.align_acf(data_1d_t,acf_1d_t,setmode) 
# See detail about modes
################################################################
##  mode       |   Ext   |         Input          |   Ext
##  -----------+---------+------------------------+---------
##  'mirror'   | 4  3  2 | 1  2  3  4  5  6  7  8 | 7  6  5
##  'reflect'  | 3  2  1 | 1  2  3  4  5  6  7  8 | 8  7  6
##  'nearest'  | 1  1  1 | 1  2  3  4  5  6  7  8 | 8  8  8
##  'constant' | 0  0  0 | 1  2  3  4  5  6  7  8 | 0  0  0
##  'wrap'     | 6  7  8 | 1  2  3  4  5  6  7  8 | 1  2  3
################################################################
def align_acf(data_1d_t, acf_1d_t, setmode):
	import numpy as np
	align_shift = convolve_1d_t_min(acf_1d_t, data_1d_t, setmode) 
	box_nbins = len(acf_1d_t[0])
	align_shift = box_nbins  - align_shift
	
	# set 0 if shifting index is at boundary
	if setmode == 'wrap':
		for index in align_shift:
			if index == box_nbins:
				index = 0 # by periodic boundary condition

	# shifting
	for iframe in range(len(data_1d_t)):
		shift_array = data_1d_t[iframe]
		data_1d_t[iframe] = np.roll(shift_array, align_shift[iframe])

	return data_1d_t
