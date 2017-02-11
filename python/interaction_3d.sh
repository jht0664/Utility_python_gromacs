#!/bin/bash
#SBATCH --job-name="TITLE"
#SBATCH --partition=pre                 # default "univ" if not specified
#SBATCH --time=24:00:00               # run time in days-hh:mm:ss
#SBATCH --nodes=1                       # require 2 nodes
#SBATCH --ntasks-per-node=1            # default 20 if this line not specified
#SBATCH --mem-per-cpu=4000              # RAM per CPU core, in MB (default 4 GB/core)

echo $(pwd)

~/anaconda3/bin/python ~/Utility/python/interaction_3d.py -b INIT -e FINAL -i traj.trr -s topol.tpr -select1 b.select -select2 hi.select -cutoff 3.0 -plot3d b-hi.plot3d.INIT > out.plot3d.INIT


