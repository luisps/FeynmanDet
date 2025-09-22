#!/bin/sh

#SBATCH --job-name=F-Det
#SBATCH --time=00:40:00
#SBATCH --partition=large-x86
#SBATCH --account=eehpc-ben-2025b09-021x
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load GCC/13.3.0

./FeynmanDet $1 $2 $3
