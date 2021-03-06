# calculation center of mass and alignment with two data sets
# input: x_t, x2_t is 1d-coordinate of atoms along time (t1, t2, ...)
#			[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
# 		 box_t is box dimension sets along time (t1, t2, ...)
#			[box_x(t0), box_x(t1), ...]
# output: com_t is center-of-mass 
#           [com_t(t0), com_t(t1), ...] 
def com_t_1d_w_data2(x_t, x2_t, box_x_t):
	import numpy as np
	print("analyze.com_t_1d:")
	# check n_frames of x_t and box_x_t
	if len(x_t) != len(box_x_t):
		raise ValueError("# of time frame is not the same for input arrarys")
	# make sure all element is within box dimension
	i_frame = 0
	for x in x_t:
		for xi in x:
			if (xi < 0) or (xi > box_x_t[i_frame]):
				xi = np.mod(xi,box_x_t[i_frame])
				print(" position is beyond box. Wrap trajectories within box dimention")
				print(" now {} within {}. ".format(xi,box_x_t[i_frame]))
		i_frame = i_frame + 1
	# make com trajectory
	com_t = np.mean(x_t,axis=1)
	com_t = np.mod(com_t,box_x_t)
	print(" COM std = {}".format(np.std(com_t)))

	# align trajectory
	i_frame = 0
	import copy
	align_x_t = copy.copy(x_t)
	align_x2_t = copy.copy(x2_t)
	for i_frame in range(len(align_x_t)):
		hello = align_x_t[i_frame] - com_t[i_frame] + (box_x_t[i_frame]/2.0)
		hello2 = align_x2_t[i_frame] - com_t[i_frame] + (box_x_t[i_frame]/2.0)
		align_x_t[i_frame] = np.mod(hello,box_x_t[i_frame])
		align_x2_t[i_frame] = np.mod(hello2,box_x_t[i_frame])
	return align_x_t, align_x2_t

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
	print("analyze.histo_t_1d_nbin:")
	import numpy as np
	# check n_frames of x_t and box_x_t
	if len(x_t) != len(box_x_t):
		raise ValueError(" # of time frame is not the same for input arrarys")
	# set number of bins (const) using input and initial box dimension
	n_bins = nbin
	# initailize variables
	n_frames = len(x_t)
	histo_t = np.zeros((n_frames, n_bins))
	bin_t = np.zeros((n_frames, n_bins+1))
	# make histogram trajectory
	i_frame = 0
	for x in x_t:
		#if (np.amax(x) > box_x_t[i_frame]) or (np.amin(x) < 0.0):
		#	print(" ##### at frame {}, index of max = {}, index of min = {}".format(i_frame,np.argmax(x),np.argmin(x)))
		#	print("Be aware of maximum value over (or less) box size, due to {} > {} (box), or {} < {} (box)##### ".format(np.amax(x),box_x_t[i_frame],np.amin(x),0.0))
		wrap_x = np.mod(x,box_x_t[i_frame])
		histo_t[i_frame], bin_t[i_frame] = np.histogram(wrap_x, bins=n_bins, range=(0,box_x_t[i_frame]))
		if np.sum(histo_t[i_frame]) != len(x):
			loss_ratio = 100*float(len(x)-np.sum(histo_t[i_frame]))/float(len(x))
			if loss_ratio > 0.1: # loss data > 0.1 %, print warning
				print(" WARNING loss of data = {} % ".format(loss_ratio))
		i_frame += 1
	return histo_t, bin_t

# we make average distribution of y values with respect to x values in a given bin_1d array
# input: x_t is 1d-coordinate of atoms along time (t1, t2, ...)
#			[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
#        y_t is properties (such as orientational parameter) of atoms along time (t1, t2, ...)
#           [[y1(t0), y2(t0), y3(t0), ...], [y1(t1), y2(t1), y3(t1)], ...]
#           where the order should be matched with x_t because order means residue id
#        bin_1d_t is bin array from np.histogram
#           
def histo_xy_t_1d_wbin(x_t, y_t, bin_1d_t):
	print("analyze.histo_xy_t_1d_wbin:")
	import numpy as np
	# check sizes
	if len(x_t) != len(y_t):
		raise ValueError(" # of time frame is not the same for input arrarys1")
	if len(y_t) != len(bin_1d_t):
		raise ValueError(" # of time frame is not the same for input arrarys2")	
	n_frames = len(x_t)
	if len(x_t[0]) != len(y_t[0]):
		raise ValueError(" # of atoms (or molecules) is not the same for input arrarys")
	n_bins = len(bin_1d_t[0]) - 1
	bin_sizes = bin_1d_t[:,1] - bin_1d_t[:,0]
	# check possible errors
	if np.amin(x_t) < 0:
		raise ValueError(" not support for negative values for x's")
	if np.amax(x_t) > np.amax(bin_1d_t[:,-1]):
		raise ValueError(" not support for values for x's beyond bin range")
	if np.amin(bin_1d_t) < 0:
		raise ValueError(" not support for negative values for bins")
	# initialize variables
	histo_t = np.zeros((n_frames, n_bins))
	# calc. histogram for y's
	for i_frame in range(n_frames):
		histo_t_count = np.zeros(n_bins,dtype=int)
		# find histogram index from x values
		x_index = (x_t[i_frame]/bin_sizes[i_frame]).astype(int)
		for i_mol in range(len(x_index)):
			if x_index[i_mol] == n_bins:
				histo_index = n_bins - 1
			else:
				histo_index = x_index[i_mol]
			histo_t[i_frame,histo_index] = histo_t[i_frame,histo_index] + y_t[i_frame,i_mol] 
			histo_t_count[histo_index] = histo_t_count[histo_index] + 1
		# average
		for ibin in range(n_bins):
			if histo_t_count[ibin] != 0:
				histo_t[i_frame,ibin] = histo_t[i_frame,ibin]/histo_t_count[ibin]
	return histo_t

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
	print("analyze.autocorr_1d_t:")
	
	import numpy as np
	from scipy import ndimage

	acf_data = np.zeros(np.shape(data_1d_t))
	itime = 0
	for data_1d in data_1d_t:
		delta_data_1d = data_1d - data_1d.mean()
		acf_data_1d = ndimage.correlate(delta_data_1d,delta_data_1d,mode=setmode)
		acf_data_1d /= (data_1d.var()*len(delta_data_1d)) # normalize
		acf_data[itime] = acf_data_1d
		itime += 1
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
	print("analyze.convolve_1d_t_min: ######### remove soon this function ###################")

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

def convolve_1d_t(data1_1d_t, data2_1d_t, setmode, minmax):
	print("analyze.convolve_1d_t:")
	import numpy as np
	from scipy import ndimage
	import copy
	# check minmax argument
	if (minmax != 'max') and (minmax != 'min'):
		raise ValueError(" Error: wrong arugment on minmax")
	# check total # element and length of datas
	data1_1d_t = np.array(data1_1d_t)
	data2_1d_t = np.array(data2_1d_t)
	# for this case, the size of data array should be the same
	if data1_1d_t.size != data2_1d_t.size:
		raise ValueError(" Error: # elements of datas are not same.")
	len_data1 = len(data1_1d_t)
	len_data2 = len(data2_1d_t)
	if len_data1 != len_data2:
		raise ValueError(" Error: length of datas are not same.")
	if len_data1*len(data1_1d_t[0]) != data1_1d_t.size \
		or len_data2*len(data2_1d_t[0]) != data2_1d_t.size:
		raise ValueError(" Error: datas are not homogeneous shape.")
	
	# Do convolution at the same time frame
	convolve_data = np.full_like(data1_1d_t,0.)
	for iframe in range(len_data1):
		convolve_temp = ndimage.convolve(data1_1d_t[iframe],data2_1d_t[iframe],mode=setmode)
		# Not neccesary normalization
		convolve_data[iframe] = copy.copy(convolve_temp)
	# keep in mind that the x-range of convolve_data starts 0 to size_data
	if minmax == 'min':
		output = np.argmin(convolve_data, axis=1)
	else:
		output = np.argmax(convolve_data, axis=1)
	output = int(len(data1_1d_t[0])/2)-output
	return output

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
		if iframe > 0:
			shift_bins = align_shift[iframe] - align_shift[iframe-1]
			if shift_bins >= 5:
				print("problem with alignment, shifting a lot by {} bins at {} iframe".format(shift_bins,iframe))
		data_1d_t[iframe] = np.roll(shift_array, align_shift[iframe]) #align_shift[0]

	return data_1d_t

def align_acf_w_data2(data_1d_t, data2_1d_t, acf_1d_t, setmode):
	import numpy as np
	print("analyze.align_acf_w_data2:")
	# acf function set all positive elements
	align_shift = convolve_1d_t(acf_1d_t, data_1d_t, setmode, 'max') 
	print(" Convolution std = {}".format(np.std(align_shift)))
	
	# shifting
	for iframe in range(len(data_1d_t)):
		shift_array = data_1d_t[iframe]
		shift_array2 = data2_1d_t[iframe]
		data_1d_t[iframe] = np.roll(shift_array, align_shift[iframe]) #align_shift[0]
		data2_1d_t[iframe] = np.roll(shift_array2, align_shift[iframe]) #align_shift[0]

	return data_1d_t, data2_1d_t

def diff_com_conv_w_data4(data11_1d_t, data12_1d_t, data21_1d_t, data22_1d_t, setmode):
	import numpy as np
	print("analyze.diff_com_conv_w_data4:")
	align_shift = convolve_1d_t(data11_1d_t,  data21_1d_t, setmode, 'max') 
	print(" diffrence shift std = {}".format(np.std(align_shift)))
	
	# shifting
	for iframe in range(len(data11_1d_t)):
		shift_array11 = data21_1d_t[iframe]
		shift_array12 = data22_1d_t[iframe]
		data21_1d_t[iframe] = np.roll(shift_array11, align_shift[iframe]) #align_shift[0]
		data22_1d_t[iframe] = np.roll(shift_array12, align_shift[iframe]) #align_shift[0]

	# difference
	diff_data1_1d_t = data11_1d_t - data21_1d_t
	diff_data2_1d_t = data12_1d_t - data22_1d_t

	return diff_data1_1d_t, diff_data2_1d_t

# fast version of einstein relation calculation using custom_func form
# input: data is a array to input 
#         [ [[x1(t1),y1(t1),z1(t1)], [x2(t1),y2(t1),z2(t1)],...], [x1(t2),y1(t2),z1(t2)], [x2(t2),y2(t2),z2(t2)],...]]
#        func1  is a defined function in python
#        mode   is (None, idx)
# output: result is a array using time delay
#         result[0] = x1(t2) - x1(t1), x2(t2) - x2(t1), ...
#         result[1] = x1(t3) - x1(t1), x2(t3) - x2(t1), ...
#     ... result[.] = x1(t3) - x1(t2), x2(t3) - x2(t2), ...
#     .. result[..] = x1(t4) - x1(t3), x2(t4) - x2(t3), ...
#       not only for x but also y and z...
# Example: def pbc(a,box): return a - box*np.around(a/box-0.5); einstein_time_relation(data,pbc,box)
def pbc(a,box):
	import numpy as np
	return a - np.around(a/box-0.5)*box

def einstein_time_relation(data,func1=None,func1_var1=None,mode=None):
	import numpy as np

	data_size = len(data)
	N = data_size*(data_size-1)//2
	idx = np.concatenate(( [0], np.arange(data_size-1,0,-1).cumsum() ))
	start, stop = idx[:-1], idx[1:]
	if len(data.shape) != 1:
		out = np.empty((N,data.shape[1]),dtype=data.dtype)
	else:
		out = np.empty(N,dtype=data.dtype)
	if func1 is None:
		for j,i in enumerate(range(data_size-1)):
			out[start[j]:stop[j]] = data[i+1:] - data[i,None]
	else:
		for j,i in enumerate(range(data_size-1)):
			out[start[j]:stop[j]] = pbc(data[i+1:] - data[i,None],func1_var1)
	if mode is None:
		return out
	elif 'idx' in mode:
		return out, idx

# generate dt list for msd, viscosity, and dynamical properties
# input:
# 	max_dt the maximum value of dt, i.e. nframes 
def gen_dt_list(max_dt):
	import numpy as np

	list_dt = []
	basic_list = np.arange(1,10)
	while True: 
		if basic_list[0] >= max_dt:
			break
		for element in basic_list:
			if element < max_dt:
				list_dt.append(element)
			else:
				break
		basic_list = basic_list*10
	
	basic_list = np.array(list_dt,dtype=int)
	return basic_list

