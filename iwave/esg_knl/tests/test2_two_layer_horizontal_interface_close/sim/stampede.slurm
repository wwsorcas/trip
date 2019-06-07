#!/bin/bash

#SBATCH -J ecsg
#SBATCH -o test.o%j
#SBATCH -e test.e%j
#SBATCH -p normal
#SBATCH -t 44:00:00
#SBATCH --mail-user=mz10@rice.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

#SBATCH -n 1
#SBATCH -N 1

# set environment variables
export OMP_NUM_THREADS=16
export OMP_PROC_BIND=close
export OMP_PLACES="{0}:16"

cd /work/02049/mz10/ECSG/tests/test2_two_layer_horizontal_interface_close/sim &&
scons

