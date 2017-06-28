import numpy as np

# set default value of tolerance of block average unless argparse does not support default
# input: tol is a float or integer
#		 deafult is default value you want to set (depends on your program)
# output: convert to a float value of tolerance, otherwise default value
# Example: args.tol = default_value(args.tol, 0.0)
def default(tol, default):
	print("blockavg.default: ")
	if tol is None:
		return float(default)
	else:
		return float(tol)

# print settings
# input: tol is a setting for printing
def print_init(tol):
	print("blockavg.print: ")
	if tol == 0.0:
		print(" Set no bloack average")
	elif tol <= 1.0: 
		print(" Tolerance for block average = %f" %tol)
	else: # for tol > 1.0
		print(" Set block length in frames = %d" %(int(tol)))
	return

def check(tol):
	print("blockavg.check: ")
	if tol < 0.0:
		raise ValueError(" Wrong input of tolerance, %f" %tol)
	elif tol > 1.0: 
		print(" The tolerance %f is assigned to the block_size, %d" %(tol, int(tol)))
		return int(tol) 
	return tol

# main program of block average of a data using tolerance of std of self or another data
# all averages of blocks should be less than (average of total data +- tolerance*std)
# (x1(t0)+x1(t1)+x1(t2)+....)/block_length = block average
# input: opt_data_t is array along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
#			block_length is block length you want
# output: new_data_t is averaged block array and block_length
# Example: number_t_1d, block_length = hjung.blockavg.main_1d(number_t_1d, unit_cells_1d, args.tol) 
# Assume: opt_data_t matrix size is homogeneous. In other words, all element array in opt_data_t is the same size
def main_1d(opt_data_1d_t, ref_data_1d, tolerance):
	print("blockavg.main_1d: ")
	tolerance = check(tolerance)

	# determine block_length
	block_length = -1
	if tolerance <= 1 and tolerance > 0:
		if ref_data_1d is None:
			print(" use itself to optimize block length")
			print(" find the max. optimal length of each x1, x2, x3, etc")
			data_t_1d = opt_data_1d_t.transpose()
			opt_lengths = []
			for data_t in data_t_1d:
				temp = opt_length_1d(data_t,tolerance) 
				opt_lengths.append(temp)
			block_length = np.amax(opt_lengths)
		else:
			block_length = opt_length_1d(ref_data_t,tolerance)
		
		if block_length == 1:
			return opt_data_1d_t, block_length
	
	elif tolerance > 1:
		block_length = int(tolerance)
	
	elif tolerance == 0:
		block_length = 1
		print(" no block averaged")
		return opt_data_1d_t, block_length
	
	else:
		raise ValueError(" tolerance error. Please check again!")
	
	if block_length <= 1:
		raise ValueError(" block length error, %d. Please check again!" %block_length)	

	# run block average with block_length
	nitems = len(opt_data_1d_t[0])
	nframes, remain = divmod(len(opt_data_1d_t),block_length)
	print(" set block length, %d" %block_length)
	print(" uses %d blocks" %nframes)
	new_data_t = np.zeros((nframes,nitems)) # new data to store average
	temp = np.zeros(nitems) # temporary data set
	i = 0
	j = 0
	for iframe in range(remain,len(opt_data_1d_t)):
		temp = np.add(temp,opt_data_1d_t[iframe])
		i += 1
		if np.mod(i,block_length) == 0:
			new_data_t[j] = temp/float(block_length)
			j += 1 
			temp = np.zeros(nitems) # reset temp
	return new_data_t, block_length

# get the optimum block length for block average using tolerance of std of data
# all averages of blocks should be less than (average of total data +- tolerance*std)
# input: data_1d is 1D-data trajectory
#		[x_t1, x_t2, x_t3...]
#			tolerance is what you want to set acceptance line
# output: opt_block_length (integer) is the optimal block length. If not found, zero.  
# Example: opt_block_length = opt_length_1d(unit_cell_1d, tolerance)
# Assume: if mod(length of data, block_length) != 0, only last trajectory is considered, not beginning.
def opt_length_1d(data_1d, tolerance):
	print("opt_length_1d: ")
	print(" To optimize, we use unit cell length on the axis you select.")
	# info before block average
	ref_avg = np.mean(data_1d)
	ref_std = np.std(data_1d)
	ref_high = ref_avg + ref_std*tolerance
	ref_low = ref_avg - ref_std*tolerance
	#print("(before block avg) Total data = %f +- %f" %(ref_avg,ref_std))
	if (100*float(ref_std)/ref_avg) < 0.05: # std is very small (<0.05%)
		print(" original std value is less than 5%. block average not necessary.\n")
		return 1
	
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
		print(" Not found optimal block length.")
		return 2
	else:
		print(" Found optimal block length, %d." %opt_block_length)
		print(" Before bloack average = %f +- %f " %(ref_avg,ref_std))
		print(" -> After bloack average = %f +- %f " %(out_avg,out_std))
		return opt_block_length


