#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on  6/5/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='NVT ensemble MC simulation with Widom-Rowlinson model')
## args
parser.add_argument('-i', '--input', default='conf.npz', nargs='?', 
	help='input .npz file')
parser.add_argument('-sizea', '--sizea', default=1.0, nargs='?', 
	help='diameter of A')
parser.add_argument('-sizeb', '--sizeb', default=1.0, nargs='?', 
	help='diameter of B')
parser.add_argument('-run', '--run', default=10.0, nargs='?', 
	help='10^12 MC steps to move particles')
parser.add_argument('-save', '--save', default=1.0, nargs='?', 
	help='10^9 MC steps to save configurations')
parser.add_argument('-o', '--output', default='confout.npz', nargs='?', 
	help='output file')
parser.add_argument('-t', '--traj', default='traj.npz', nargs='?', 
	help='traj file')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# default for args
args.input = args.input if args.input is not None else 'conf.npz'
args.sizea = args.sizea if args.sizea is not None else 1.0
args.sizea = float(args.sizea)
args.sizeb = args.sizeb if args.sizeb is not None else 1.0
args.sizeb = float(args.sizeb)
args.step_run = args.step_run if args.step_run is not None else 10.0
args.step_run = float(args.run*(10**12))
args.step_save = args.step_save if args.step_save is not None else 1.0
args.step_save = float(args.run*(10**9))
args.output = args.output if args.output is not None else 'confout.npz'
args.traj = args.traj if args.traj is not None else 'traj.npz'

# numpy double precision
import numpy as np
args.sizea = np.float_(args.sizea)
args.sizeb = np.float_(args.sizeb)
args.step_run = np.float_(args.step_run)
args.step_save = np.float_(args.step_save)

## timer
import time
start_clock = time.clock() # process time
start_wall = time.time() # wall time

# load input file
print("="*30)
print("load input file, %s" %args.input)
inputfile = np.load(args.input) # str = nmola, nmolb, coord, box
nmola = inputfile['nmola']
nmolb = inputfile['nmolb']
ntot = nmola + nmolb
coordinates = inputfile['coord']
box = inputfile['box']
inputfile.close()

# check input file
print("Assume coordinates has A molecules, then B molecules in order")
print("Total number of molecules = %s" %ntot)
if len(coordinates) != ntot:
	raise OSError("size of coordinates does not match total # molecules.")
print("Box = %s" %(box))

# mapping cells to reduce distance calculation
cell_size = max(args.sizea,args.sizeb)*1.2
ncell = np.trunc((box/cell_size))
print("cell matrix = %s" %ncell)
cell_xyz = box/ncell
ntot_cell = np.prod(ncell)
cell_map = np.zeros((ntot_cell,27),dtype=np.int)
for ic in range(ntot_cell):
	icx = ic%ncell[0]
	icy = ((ic-icx)/ncell[0])%ncell[1]
	icz = (ic-icx-ncell[0]*icy)/ncell[1]/ncell[0]
	id = 0
	for ix in range(icx-1,3):
		if ix = -1:
			ix = ncell[0]-1
		if ix >= ncell[0]:
			ix = ix - ncell[0]
		for iy in range(icy-1,3):
			if iy = -1:
				iy = ncell[1]-1
			if iy >= ncell[1]:
				iy = iy - ncell[1]
			for iz in range(icz-1,3):
				if iz = -1:
					iz = ncell[2]-1
				if iz >= ncell[2]:
					iz = iz - ncell[2]

				ict = ncell[0]*ncell[1]*iz +ncell[0]*iy+ix
				cell_map[ic][id] = ict
				id = id + 1

# mapping molecuels
mol_map = np.zeros(ntot,dtype=np.int)
for i in range(ntot):
	icell = trunc(coordinates[i][0]/cell_xyz[0]) + 
	trunc(coordinates[i][1]/cell_xyz[1])*ncell[0] + 
	trunc(coordinates[i][2]/cell_xyz[2])*ncell[0]*ncell[1]
	mol_map[i] = icell

# run mc simulations
nruns = 0
naccept = 0
# pick a type of molecules
mol_try = np.randint(0,ntot)

# return a cell number of the coord using # cells on x,y,z and length
def get_icell(coord,ncell,cell_xyz):
	return trunc(coord[0]/cell_xyz[0]) + 
	trunc(coord[1]/cell_xyz[1])*ncell[0] + 
	trunc(coord[2]/cell_xyz[2])*ncell[0]*ncell[1]

new_coord = [np.random.random_sample()*box[0],
		np.random.random_sample()*box[1], 
		np.random.random_sample()*box[2]]
icell = get_icell(new_coord,ncell,cell_xyz)
neighbor_cell = cell_map[icell]
list_mol = []
for i in range(ntot):
	for i_neighbors in neighbor_cell:
		if mol_map[i] == i_neighbors:
			list_mol.append(i)
# check overlap
for i in range(nmola):
		dxyz = coordinates[i] - coordinates[new_index]
		dxyz = dxyz - box * np.around(dxyz/box)
		if np.linalg.norm(dxyz) < distance:
			return 1 # overlap!
	return 0 # success for insertion

for i in range(args.nmola):
	
curr_index = args.nmola
ntry = 1
nmod = 10
while (curr_index < ntot):
	if ntry%nmod == 0:
		print("%d th trials (%s/%s)th particle" %(ntry,curr_index,ntot))
		if (ntry/nmod)%nmod == 0:
			nmod = nmod*nmod
	if ntry > args.maxtry:
		print("Hard to insert because ntry > maxtry.")
		print("I made initial coordinates with %d out of %d molecules" %(curr_index-1,ntot))
		break
	coordinates[curr_index] = [np.random.random_sample()*box[0],
		np.random.random_sample()*box[1], 
		np.random.random_sample()*box[2]]
	dist = 0.5*(args.sizea + args.sizeb)
	success = overlap(curr_index,coordinates,args.nmola,dist,box)
	if success == 0:
		curr_index = curr_index + 1
	ntry = ntry + 1

# save initial coordinates
print("="*30)
print("Saving OutputFile...")
output_file = open(args.output, 'w')
output_file.write('# generated by initial.py')
output_file.write('%d' %(ntot))
for i in range(min(curr_index,ntot)):
	if i < args.nmola:
		output_file.write('{0:5d}{1:<5s}{2:5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}'.format(
			i+1,'LJA','A',i+1,coordinates[i][0],coordinates[i][1],coordinates[i][2]))
	else:
		output_file.write('{0:5d}{1:<5s}{2:5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}'.format(
			i+1,'LJB','B',i+1,coordinates[i][0],coordinates[i][1],coordinates[i][2]))
output_file.write('{0:10.5f}{1:10.5f}{2:10.5f}'.format(box[0],box[1],box[2]))
output_file.close()

print("="*30)
print(time.clock() - start_clock, "seconds process time")
print(time.time() - start_wall, "seconds wall time")

#class MaxError(Error):
#	def __init__(self, expression, message):
#    	self.expression = expression
#    	self.message = message
#
## assign particles until no overlapping
#def overlap(maxtry,new_index,coordinates,distance,unit_cells):
#	ntry = 0
#	try:
#		return overlap_recur(ntry,maxtry,new_index,coordinates,distance,unit_cells)
#	except MaxError:
#		print ('You need to use alternataive way for insertion to avoid overlapping.')
#	else:
#		print('Different error in the function, overlap.')
#	return
#
#def overlap_recur(ntry,maxtry,new_index,coordinates,distance,unit_cells):
#	if ntry > maxtry:
#		raise MaxError('# trials, %d > Maxtry, %d!!' %(ntry,maxtry))
#		return
#	import numpy as np
#    for i in range(len(coordinates)):
#		if np.linalg.norm(coordinates[i]-coordinates[new_index]) < distance:
#			print("Overlap!")


