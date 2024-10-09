#!/bin/bash -l
#SBATCH --job-name=paths
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/paths_%A.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/paths_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=12:00:00

execdir=/cellar/users/snwright/Git/rare_common/carva
srun -l python $execdir/find_shortest_paths.py


