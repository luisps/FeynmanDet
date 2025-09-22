#!/bin/sh

sbatch <<EOT
#!/bin/sh

#SBATCH --job-name=Feynman-C$1-A$3-T$4
#SBATCH --output=Feynman-C$1-A$3-T$4-out.txt
#SBATCH --error=Feynman-C$1-A$3-T$4-error.txt
#SBATCH --time=06:00:00
#SBATCH --partition=normal-x86
#SBATCH --account=I20250002X
#SBATCH --nodes=1
#SBATCH --ntasks=$4

module load GCC/13.3.0

./FeynmanDet-omp $1 $2 $3 $4
EOT
