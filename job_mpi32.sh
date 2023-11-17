#!/bin/bash

#SBATCH --job-name=MyJob
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt

#SBATCH --partition=cpu_test
#SBATCH --account=ams301

#SBATCH --ntasks=32
#SBATCH --time=00:10:00

## load modules
module load cmake
module load gcc
module load gmsh
module load openmpi

## execution
mpirun -display-map ${SLURM_SUBMIT_DIR}/solver

