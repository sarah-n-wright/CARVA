#!/bin/bash -l
#SBATCH --job-name=matrix
#SBATCH --output /cellar/users/snwright/Data/RareCommon/matrix_%A.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/matrix_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=01:00:00

OUTDIR=$1

srun -l python get_heat_matrix.py $1



