import numpy as np

# global variable random set
rndset = np.random.RandomState(1985)

# print process percentage
# input: istep is the current nstep
#		 nstep is the final nstep
#        modn is mode number to print
# output: modn is mode number to print
# Example: modn = print_process(istep, tstep, modn)
def print_process(text, istep, nstep, modn):
	#print('{0} {1} {2}'.format(istep, nstep, modn))
	if istep%modn == 0:
		print("... {0}/{1} th ({2:.0%}) for {3}...".format(istep,nstep,istep/nstep,text))
		if (istep/modn)%10 == 0:
			modn = modn*10
	return modn

# set default value of tolerance of block average unless argparse does not support default
# input: nmon is # total monomers (or beads).
#		 ndens is number density of monomers (or beads)
#        ratio is (x, y, z) ratio (= [int, int, int])
# output: initial box size array, [box_x, box_y, box_z]
# Example: lj_box = lj_box_init_w_ndens_ratio(nmon, ndens, [1, 1, 5])
def lj_box_init_w_ndens_ratio(nmon, ndens, ratio):
	print("lj_box_init_w_ndens_ratio: ")
	print(" total number of LJ particles = %s" %nmon)
	boxl = (float(nmon)/ndens/np.prod(ratio))**(1.0/3.0) 
	box = np.multiply(ratio,boxl)
	print(" box array = [{0:8.3f}, {1:8.3f}, {2:8.3f}]".format(box[0],box[1],box[2]))
	return box

# get a pre-separated coordinate using lattice insertion
# To reduce cost to check overlapping, use cell list with 0 (empty), 1 (occupied A), and 2 (oocupied B).
# input: nmol_a is # A chins
#		 mon_a is # monomers of a A chain
#        cell_length is the length of lattice (cubic) cell
#        box is box dimension [box_x, box_y, box_z]
#        frac_a is fraction of A in one of coexisting two phases
#        seed is seed number for random number generator
# output: coordinates of beads ((x1,y1,z1),(x2,y2,z2),...)
# Example: coordinates = hjung.initial.insertion_lattice_sep(args.nmola, args.mona, args.nmolb, args.monb, cell_length, box, args.frac, args.seed)
def insertion_lattice_sep(nmol_a, mon_a, nmol_b, mon_b, guess_cell_length, box, frac_a, seed):
	print("insertion_lattice_sep:")
	
	# assign random set 
	global rndset
	rndset = np.random.RandomState(seed)

	if box.shape != (3,):
		raise ValueError(" the input 'box' is not (3,) array!")
	# check validity of lattice insertion method
	cell = np.floor(np.divide(box,guess_cell_length)).astype(int)
	print(" cell dimension = [{0}, {1}, {2}]".format(cell[0],cell[1],cell[2]))
	cell_length = np.divide(box,cell.astype(float))
	print(" cell length = [{0:8.3f}, {1:8.3f}, {2:8.3f}]".format(cell_length[0],cell_length[1],cell_length[2]))
	ncell = np.prod(cell)
	nmtot = nmol_a*mon_a + nmol_b*mon_b
	nctot = nmol_a + nmol_b
	if ncell < nmtot:
		raise ValueError(" Not possible to use lattice insertion because #particles > #cells. Reduce cell_length.")

	# pre-determine number of A and B chains in a phase
	sep_maxa = int(nmol_a*frac_a)
	print(" one phase has {0} of A chains".format(sep_maxa))
	sep_maxb = int(nmol_b*(1.0-frac_a))
	print(" also, the phase has {0} of B chains".format(sep_maxb))
	
	# lattice insertion
	cell_list = np.zeros(cell,dtype=int) # Cell list with 0 (empty), 1 (occupied A), and 2 (oocupied B).
	modn = 10
	icell_a, cell_list = insert_polymer_sep("monomers of success insertion A", 1, sep_maxa, nmol_a, mon_a, box, cell_list, modn)
	icell_b, cell_list = insert_polymer_sep("monomers of success insertion B", 2, sep_maxb, nmol_b, mon_b, box, cell_list, modn)

	# save coordinate
	coordinates = np.zeros((nmtot, 3))
	index_coord = 0
	for itr_icell_a in icell_a:
		coordinates[index_coord] = (itr_icell_a+0.50)*cell_length
		index_coord += 1
	for itr_icell_b in icell_b:
		coordinates[index_coord] = (itr_icell_b+0.50)*cell_length
		index_coord += 1
	if index_coord != nmtot:
		raise RuntimeError(" not match index of total #particles.")

	return coordinates

# insert polymer chains of a component
# input: text_a is for print_process
#		 monomer_type is numbering of monomer type, int.
#		 sep_max is # monomers in a pre-separated phase 
#        nmol is # chains of the polymer
#		 mon is # monomers of the polymer
#        box is box dimension [box_x, box_y, box_z]
#        cell_list is the 3d array with information for empty or occupied by monomer_type
# output: coordinates of beads ((x1,y1,z1),(x2,y2,z2),...) 
#         cell_list is the 3d array with information for empty or occupied by monomer_type
# Example: coordinates, cell_list = insert_component("Trials of Insertion A", sep_maxa, nmol_a, box, cell_length, cell, cell_list, coordinates, 10, seed)
def insert_polymer_sep(text, monomer_type, sep_max, nmol, mon, box, cell_list, moden):
	print("insert_component:")
	polymer = np.zeros((nmol*mon,3),dtype=int)
	sep_layer1 = np.array([1.0, 1.0, 0.50]) # for ratio of phases
	sep_layer2 = np.array([0.0, 0.0, 0.50]) # for translation
	imol = 0
	curr_index = 0
	while (imol < nmol):
		# print process
		moden = print_process(text,imol,nmol,moden)
		generate_mol = False
		
		while generate_mol is False:
		# make pre-separated phase, left and right
			if imol < sep_max:
				monomer_seed = generate_seed_sep('lattice',sep_layer1,np.array([0.0, 0.0, 0.0]),cell_list.shape)
			else:
				monomer_seed = generate_seed_sep('lattice',sep_layer1,sep_layer2,cell_list.shape)
				
			# make a cell list and check validity of generated a cell list
			#monomer_seed = np.trunc(monomer_seed*box/cell_length) 
			if check_monomer(monomer_seed, cell_list):
				cell_list = update_monomer(monomer_seed, monomer_type, cell_list)
			else:
				continue
		
			# start growth of polymer
			success, icell_new_chain, new_cell_list = insert_chain(monomer_seed, monomer_type, mon, cell_list) # 1 = monomer type A
			if success is True:
				# update cell list
				cell_list = new_cell_list
				# save the coordinate of the new chain
				for itr_new_chain in icell_new_chain:
					polymer[curr_index] = itr_new_chain
					curr_index += 1
				imol += 1
				generate_mol = True
			else:
				# get a new monomer seed again
				#print(" fail to grow. redo generate seed.")
				continue

	# check index 
	if curr_index != nmol*mon:
		raise ValueError(" total # particles of the component is not the same with your setting and algorithm.")

	return polymer, cell_list

# generator a monomer seed for growth of polymer 
# input: sep_layer_ratio is separated layer ratio (all elements are <= 1), 
#		 sep_layer_trans is separated layer translation (all elements are <= 1),
#        max_val is maximum value (exclusive) of random sample for lattice mode, [max_x, max_y, max_z] for lattice mode
# output: monomer_seed is random position of [0,1) with separation presets
# Example: success, icell_new_chain, new_cell_list = insert_chain(iclist, 1, mon_a, cell_list, seed)
#################################################################################
# For example, ratio [1.0, 1.0, 0.50] and trans [0.0, 0.0, 0.50]
# no separation on x and y, but left 50% of box space is picked for the monomer                  
# no translation of phases on x and y, but the left 50% of box_z is translated by 0.50 amount of 50% length of box_z
# In other words, the monomer seed is picked in the space of right 50% of box_z.
###################################################################################
def generate_seed_sep(mode,sep_layer_ratio,sep_layer_trans,max_val):
	#print("generate_seed_sep:")
	global rndset
	monomer_seed = np.zeros(3)

	# get random monomer_seed depending on separation presets
	monomer_seed = rndset.random_sample(3)
	monomer_seed = monomer_seed * sep_layer_ratio
	monomer_seed = monomer_seed + sep_layer_trans
	if mode == 'lattice':
		if max_val is None:
			raise ValueError(" need max_val argument.")
		monomer_seed = np.array(monomer_seed*max_val, dtype=int)
	else:
		raise ValueError(" wrong mode argument.")

	return monomer_seed

# check validity of monomer
# input: monomer is a given position in cell_list
#        cell_list is the 3d array with information for empty or occupied by monomer_type
# output: true if no overlapping. False if already occupied
# Example: success = check_monomer(monomer,cell_list)
def check_monomer(monomer,cell_list):
	#print("check_seed:")
	# check overlapping of monomer
	if cell_list[monomer[0]][monomer[1]][monomer[2]] != 0:
		# already occupied the monomer position
		return False
	# no overlaps
	return True

# update monomer in cell_list
# input: monomer is a given position in cell_list
#        monomer_type is type of monomer (int.)
#        cell_list is the 3d array with information for empty or occupied by monomer_type
# output: cell_list is a updated cell_list
# Example: cell_list = update_monomer(monomer,monomer_type,cell_list)
def update_monomer(monomer,monomer_type,cell_list):
	cell_list[monomer[0]][monomer[1]][monomer[2]] = monomer_type
	return cell_list

# grow a chain 
# input: monomer_seed is icell array [int, int, int] of first monomer of the chain
#		 monomer_type is numbering of monomer type, int.
#        chain_length is degree of polymerization of the chain
#        cell_list is the 3d array with information for empty or occupied by monomer_type
# output: cell_list (updated) 
#         new_chain is the icell array for all monomers of the chain 
#         True if success. Otherwise, False
# Example: success, icell_new_chain, new_cell_list = insert_chain(iclist, 1, mon_a, cell_list, seed)
def insert_chain(monomer_seed,monomer_type,chain_length,cell_list):
	#print("insert_chain:")
	global rndset
	new_chain = np.zeros((chain_length,3),dtype=int)
	new_chain[0] = monomer_seed
	#print(" monomer seed = {0}".format(monomer_seed))

	ibead = 1
	curr_monomer = monomer_seed
	max_trials = 1000*chain_length
	ntrials = 0
	while (ibead < chain_length):
		# check # trials
		ntrials = ntrials + 1
		if ntrials > max_trials:
			break
		# determine direction
		rnd_direction = rndset.random_sample(1)*6
		if rnd_direction < 1.0:
			direc = np.array([1, 0, 0],dtype=int)
		elif rnd_direction < 2.0:
			direc = np.array([-1, 0, 0],dtype=int)
		elif rnd_direction < 3.0:
			direc = np.array([0, 1, 0],dtype=int)
		elif rnd_direction < 4.0:
			direc = np.array([0, -1, 0],dtype=int)
		elif rnd_direction < 5.0:
			direc = np.array([0, 0, 1],dtype=int)
		elif rnd_direction < 6.0:
			direc = np.array([0, 0, -1],dtype=int)
		# check overlapping
		try_index = curr_monomer + direc
		try_index = np.mod(try_index,cell_list.shape) # periodic boundary condition
		
		if check_monomer(try_index,cell_list):
			#print(' before {0}'.format(cell_list[try_index[0]][try_index[1]][try_index[2]]))
			#print(' next monomer = {0}'.format(try_index))
			cell_list = update_monomer(try_index,monomer_type,cell_list)
			#print(' after {0}'.format(cell_list[try_index[0]][try_index[1]][try_index[2]]))
			new_chain[ibead] = try_index
			curr_monomer = try_index
			ibead = ibead + 1
		else:
			continue
		
	if ntrials > max_trials:
		#print(" fails to grow.")
		return False, None, None
	else:
		# check overlapping within itself
		for imonomer in range(new_chain.shape[0]):
			for jmonomer in range(imonomer+1,new_chain.shape[0]):
				if np.all(np.equal(new_chain[imonomer],new_chain[jmonomer])):
					print(new_chain[imonomer])
					print(new_chain[jmonomer])
					raise RuntimeError(" algorithm problem due to overlapping!")

		#print(" success!")
		return True, new_chain, cell_list
