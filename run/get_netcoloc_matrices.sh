#!/bin/bash -l
#SBATCH --job-name=matrix
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/matrix_%A.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/matrix_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=01:00:00

OUTDIR=/cellar/users/snwright/Data/RareCommon/inputs/
execdir=/cellar/users/snwright/Git/rare_common/carva
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
name=pcnet2_0
srun -l python $execdir/get_heat_matrix.py $OUTDIR $uuid $name



