#!/bin/bash

#SBATCH -N 1
#SBATCH -p compute
#SBATCH -t 8:00:00 
#SBATCH -U mschpc
#SBATCH -J diffusion

# source the module commands
source /etc/profile.d/modules.sh

# load the modules used to build the xhpl binary
module load cports gcc/4.9.3-gnu openmpi/1.8.6-gnu4.9.3 

# run it
##
for size in 960 1920 2880 3840 4800
do
    # program outputs one line : " nproc   time  "
    mpirun -n 7 par_refine -r $size -i 10000 -q 2 -p 2 >> "refined_timing${size}.dat" 
done
