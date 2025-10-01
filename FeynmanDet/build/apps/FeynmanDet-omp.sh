#!/bin/sh

#SBATCH --job-name=F-omp
#SBATCH --time=12:10:00
#SBATCH --partition=large-x86
#SBATCH --account=eehpc-ben-2025b09-021x
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load GCC/13.3.0
export OMP_PROC_BIND=true
export OMP_PLACES=cores

./FeynmanDet-omp $1 $2 $3 $4
