#!/bin/bash

#SBATCH -N 8        # 8 cores = 1 node on lonsdale
#SBATCH -p compute
#SBATCH -t 0:6:30 
#SBATCH -U mschpc
#SBATCH -J diffusion

# source the module commands
source /etc/profile.d/modules.sh

# load the modules used to build the xhpl binary
module load cports gcc/4.9.3-gnu openmpi/1.8.6-gnu4.9.3 

# run it
mpirun -n 63 par_refine -r 4800 -i 10000 -q 6 -p 6 -n 57
