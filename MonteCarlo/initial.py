#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on  6/1/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='Initial coordiante of particles for MC simulation')
## args
#parser.add_argument('-i', '--input', default='init.ic', nargs='?', 
#	help='input file')
parser.add_argument('-na', '--nmola', nargs='?', type=int,
	help='# particles of A')
parser.add_argument('-nb', '--nmolb', nargs='?', type=int,
	help='# particles of B')
parser.add_argument('-r', '--ratio', default=5.0, nargs='?', type=float,
	help='ratio of box-z/box-x (box-x = box-y)')
parser.add_argument('-sep', '--sep', default='NO', nargs='?', type=str,
        help='pre-separation YES/NO')
parser.add_argument('-fr', '--frac', default=1.0, nargs='?', type=float,
	help='number fraction of A of one phase if -sep YES (')
parser.add_argument('-d', '--dens', nargs='?', type=float,
	help='number density')
parser.add_argument('-sa', '--sizea', default=1.0, nargs='?', type=float,
	help='diameter of A')
parser.add_argument('-sb', '--sizeb', default=1.0, nargs='?', type=float,
	help='diameter of B')
parser.add_argument('-mt', '--maxtry', default=0, nargs='?', type=int,
	help='attemps for random insertion (if zero, do lattice insertion)')
parser.add_argument('-fm', '--format', default='MC', nargs='?', type=str,
	help='Save in fortran MC format (MC), .npz format (NPZ), or .gro format (GRO)')
parser.add_argument('-o', '--output', default='init', nargs='?', type=str,
	help='output file (exclude extension name) ')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# default for args
#args.input = args.input if args.input is not None else 'init.ic'
if args.sep != 'YES' and args.sep != 'NO':
	raise ValueError("Wrong argument for pre-separation option")
if args.sep == 'NO' and args.frac != 1.0:
	raise ValueError("-sep and -fr not matched or not set")
if args.sep == 'YES':
	if args.nmola != args.nmolb:
		raise ValueError("Not support the different #particles of A and B")
if args.format != 'MC' and args.format != 'NPZ' and args.format != 'GRO':
	raise ValueError("Wrong argument for format!")
if args.format == 'MC':
	args.output = args.output+'.ic'
elif args.format == 'NPZ':
	args.output = args.output+'.npz'
elif args.format == 'GRO':
	args.output = args.output+'.gro'

# numpy double precision
import numpy as np
import sys
sys.path.append('/home/htjung/Utility/python/')
import hjung
from hjung import *

## timer
start_proc, start_prof = hjung.time.init()

# determine box size
print("="*30)
ntot = args.nmola + args.nmolb
print("Total number of molecules = %s" %ntot)
boxl = (float(ntot)/args.dens/args.ratio)**(1.0/3.0) 
box = np.array((boxl, boxl, boxl*args.ratio))
print("Box = %s" %(box))

# attemp to insertion only for Widom-Rowlinson Mixture
coordinates = np.zeros((ntot, 3))

def overlap(new_index,coordinates,nmola,distance,box):
	import numpy as np
	for i in range(nmola):
		dxyz = coordinates[i] - coordinates[new_index]
		dxyz = dxyz - box * np.floor(dxyz/box)
		if np.linalg.norm(dxyz) < distance:
			return 1 # overlap!
	return 0 # success for insertion

print("="*30)
print("Start Insertion")
maxa = int(args.nmola*args.frac)
maxb = int(args.nmolb*args.frac)
if args.maxtry > 0:
	# try random insertion
	for i in range(args.nmola):
		if args.sep == 'YES' and i < maxa:
			# if you set specific fraction and pre-separation
			coordinates[i] = [np.random.random_sample()*box[0],
				np.random.random_sample()*box[1], 
				np.random.random_sample()*0.50*box[2]]
		elif args.sep == 'YES' and i >= maxa:
			coordinates[i] = [np.random.random_sample()*box[0],
				np.random.random_sample()*box[1], 
				(np.random.random_sample()*0.50+0.50)*box[2]]
		else:
			# if you set random
			coordinates[i] = [np.random.random_sample()*box[0],
				np.random.random_sample()*box[1], 
				np.random.random_sample()*box[2]]
		
	curr_index = args.nmola
	ntry = 1
	nmod = 10
	while (curr_index < ntot):
		if ntry%nmod == 0:
			print("%d th random trials (%s/%s)th particle" %(ntry,curr_index,ntot))
			if (ntry/nmod)%10 == 0:
				nmod = nmod*10
		if ntry > args.maxtry:
			print("Hard to insert because ntry > maxtry.")
			print("I made initial coordinates with %d out of %d molecules" %(curr_index-1,ntot))
			break
		if args.sep == 'YES' and curr_index < args.nmola+maxb:
			coordinates[curr_index] = [np.random.random_sample()*box[0],
				np.random.random_sample()*box[1], 
				(np.random.random_sample()*0.50+0.50)*box[2]]
		elif args.sep == 'YES' and curr_index >= args.nmola+maxb:
			coordinates[curr_index] = [np.random.random_sample()*box[0],
				np.random.random_sample()*box[1], 
				np.random.random_sample()*0.50*box[2]]
		else:
			coordinates[curr_index] = [np.random.random_sample()*box[0],
				np.random.random_sample()*box[1], 
				np.random.random_sample()*box[2]]
		dist = 0.5*(args.sizea + args.sizeb)
		success = overlap(curr_index,coordinates,args.nmola,dist,box)
		if success == 0:
			curr_index = curr_index + 1
		ntry = ntry + 1
else:
	# try lattice insertion
	maxsize = max(args.sizea,args.sizeb)
	ncellx = np.int(np.floor(box[0]/maxsize))
	ncelly = np.int(np.floor(box[1]/maxsize))
	ncellz = np.int(np.floor(box[2]/maxsize))
	ncell = ncellx*ncelly*ncellz

	if ncell < ntot:
		raise ValueError("Not possible to use lattice insertion because #particles > #cells")

	occupy_cell = np.zeros((ncellx,ncelly,ncellz),dtype=int)
	i = 0
	ntry = 1
	nmod = 10
	print("Try Insertion of A")
	while (i < args.nmola):
		if ntry%nmod == 0:
			print("%d th lattice trials (%s/%s)th particle" %(ntry,i,ntot))
			if (ntry/nmod)%nmod == 0:
				nmod = nmod*nmod
		
		icx = np.trunc(np.random.random_sample()*box[0]/maxsize)
		icy = np.trunc(np.random.random_sample()*box[1]/maxsize)
		if args.sep == 'YES' and i < maxa:
			icz = np.trunc(np.random.random_sample()*0.50*box[2]/maxsize)
		elif args.sep == 'YES' and i >= maxa:
			icz = np.trunc((np.random.random_sample()*0.50+0.50)*box[2]/maxsize)
		else:
			icz = np.trunc(np.random.random_sample()*box[2]/maxsize)
		if icx < ncellx and icy < ncelly and icz < ncellz:
			randx = (icx+0.5)*maxsize
			randy = (icy+0.5)*maxsize
			randz = (icz+0.5)*maxsize
			coordinates[i] = [randx,randy,randz]
			occupy_cell[np.int(icx)][np.int(icy)][np.int(icz)] = 1 
			i = i + 1
		ntry = ntry + 1

	curr_index = args.nmola
	ntry = 1
	nmod = 10
	print("Try Insertion of B")
	while (curr_index < ntot):
		if ntry%nmod == 0:
			print("%d th lattice trials (%s/%s)th particle" %(ntry,curr_index,ntot))
			if (ntry/nmod)%nmod == 0:
				nmod = nmod*nmod
		icx = np.trunc(np.random.random_sample()*box[0]/maxsize)
		icy = np.trunc(np.random.random_sample()*box[1]/maxsize)
		if args.sep == 'YES' and curr_index < args.nmola+maxb:
			icz = np.trunc((np.random.random_sample()*0.50+0.50)*box[2]/maxsize)
		elif args.sep == 'YES' and curr_index >= args.nmola+maxb:
			icz = np.trunc(np.random.random_sample()*0.50*box[2]/maxsize)
		else:
			icz = np.trunc(np.random.random_sample()*box[2]/maxsize)
		randx = (icx+0.5)*maxsize
		randy = (icy+0.5)*maxsize
		randz = (icz+0.5)*maxsize
		coordinates[curr_index] = [randx,randy,randz]
		ntry = ntry + 1
		if icx >= ncellx or icy >= ncelly or icz >= ncellz:
			continue
		elif occupy_cell[np.int(icx)][np.int(icy)][np.int(icz)] == 0:
			curr_index = curr_index + 1

# save initial coordinates
print("="*30)
print("Saving OutputFile...")
if args.format == 'NPZ':
	# array_argument = nmola, nmolb, coord, box
	np.savez(args.output,nmola=args.nmola,nmolb=args.nmolb,coord=coordinates,box=box)
elif args.format == 'GRO':
	# gromacs version
	output_file = open(args.output, 'w')
	output_file.write('# generated by initial.py\n')
	output_file.write('%d\n' %(ntot))
	for i in range(min(curr_index,ntot)):
		if i < args.nmola:
			output_file.write('{0:5d}{1:<5s}{2:5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n'.format(
				i+1,'LJA','A',i+1,coordinates[i][0],coordinates[i][1],coordinates[i][2]))
		else:
			output_file.write('{0:5d}{1:<5s}{2:5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n'.format(
				i+1,'LJB','B',i+1,coordinates[i][0],coordinates[i][1],coordinates[i][2]))
	output_file.write('{0:10.5f}{1:10.5f}{2:10.5f}\n'.format(box[0],box[1],box[2]))
	output_file.close()
elif args.format == 'MC':
	# fortran MC version 'composite.ic'
	output_file = open(args.output, 'w')
	output_file.write('{}  #NUMBER OF B particle \n'.format(args.nmolb))
	output_file.write('{}  #NUMBER OF A particle \n'.format(args.nmola))
	output_file.write('{}  #SIZE OF B \n'.format(args.sizeb))
	output_file.write('{}  #SIZE OF A \n'.format(args.sizea))
	output_file.write('{} {} {}  # BOX SIZE \n'.format(box[0],box[1],box[2]))
	# coordinates of A, then B
	for i in range(ntot):
		output_file.write('{} {} {} {} {}\n'.format(i,i,coordinates[i][0],coordinates[i][1],coordinates[i][2]))
	output_file.close()
else:
	raise RuntimeError("Sometime wrong!")

## timer
hjung.time.end_print(start_proc, start_prof)
