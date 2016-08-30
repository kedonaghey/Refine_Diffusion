#!/bin/bash

#SBATCH -N 2        # 8 cores = 1 node on lonsdale
#SBATCH -p debug
#SBATCH -t 0:30:00 
#SBATCH -U mschpc
#SBATCH -J diffusion

# source the module commands
source /etc/profile.d/modules.sh

# load the modules used to build the xhpl binary
module load cports gcc/4.9.3-gnu openmpi/1.8.6-gnu4.9.3 

# run it
mpirun -n 14 par_refine -r 20 -i 3 -q 2 -p 4
