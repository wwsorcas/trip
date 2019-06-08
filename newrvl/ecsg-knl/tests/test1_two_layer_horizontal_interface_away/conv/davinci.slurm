#!/bin/bash
#SBATCH --job-name=ecsg
#SBATCH --account=trip
#SBATCH --partition=trip
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=70:00:00
#SBATCH --mail-user=mz10@rice.edu
#SBATCH --mail-type=ALL


echo "My job ran on:"
echo $SLURM_NODELIST

export OMP_NUM_THREADS=16
export OMP_PROC_BIND=TRUE
export OMP_PLACES="{0}:16"

cd $SHARED_SCRATCH/$USER/ecsg/tests/test1_two_layer_horizontal_interface_away/conv
scons

