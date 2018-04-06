## polymer section

# check connectivity of atoms of a polymer
#   sometimes, "a broekn polymer" may happen due to periodic boundary condition
# input:
#   mda_system mda.Universe(structure,input)
#   select selection command
#   nmol   number of polymer molecules
#   cutoff cut-off distance to check connectivity
#   mode   check only the first frame unless 'all'
def check_traj_connectivity(mda_system,select,nmol,cutoff,mode):
	print("polymer.check_traj_connectivity:")
	import MDAnalysis as mda
	import random
	from scipy.spatial.distance import euclidean
	u = mda_system
	n_frames = len(u.trajectory)
	select_mol = u.select_atoms(select)
	#print(len(select_mol),nmol)
	if len(select_mol)%nmol != 0:
		raise ValueError(" wrong # molecules, (args.nmol, select_mol) {} {} ".format(nmol, len(select_mol)))
	n_deg = int(len(select_mol)/nmol)
	print(" assume all molecules has {} atoms".format(n_deg))
	dist_min = 100.0
	dist_max = 0
	# read trajectory
	if 'all' in mode:
		set_ts = u.trajectory
		print(" active ALL mode")
	elif 'random' in mode:
		pick_frame = random.randrange(0,n_frames)
		set_tx = u.trajectory[pick_frame]
		print(" active RANDOM mode")
	else:
		set_tx = u.trajectory[0]
		print(" active FIRST mode")
	# check
	i_frame = 0
	for ts in set_tx:
		for i_mol in range(nmol):
			# check the validity of atomic positions 
			for i_atom in range(n_deg-1):
				inum = i_mol*n_deg+i_atom
				dist = euclidean(select_mol.positions[inum], select_mol.positions[inum+1])
				if dist > cutoff:
					print(" maybe due to the wrapped trajectory setting or too small cutoff.")
					print(" please check wrapping option like gmx trjconv -pbc mol for Gromacs trajectory.")
					raise RuntimeError("[{}th frame] {}th polymer ({}th atom) dist. = {} > cutoff {}".format(i_frame,i_mol,i_atom, dist, cutoff))
				if dist > dist_max:
					dist_max = dist
				if dist < dist_min:
					dist_min = dist
		i_frame = i_frame + 1
	print(" passed! with distance [{:.3f},{:.3f}] under {:.3f} cut-off".format(dist_min,dist_max,cutoff))
	return
