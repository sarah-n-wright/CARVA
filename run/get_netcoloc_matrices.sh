#!/bin/bash -l
#SBATCH --job-name=matrix
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/matrix_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/matrix_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=01:00:00
#SBATCH --array=0

OUTDIR=/cellar/users/snwright/Data/RareCommon/inputs/
execdir=/cellar/users/snwright/Git/rare_common/carva
#uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
#name=pcnet2_0
#uuid='8b4b54fa-e87d-11ee-9621-005056ae23aa'
#name=pcnet2_2
#uuid='db921c7d-f144-11ee-9621-005056ae23aa'
#name=humannet
uuid='5f5da339-f14a-11ee-9621-005056ae23aa'
name=stringHC
filter=Score_700

#uuid_list=$1
#net_names=$2
#uuids=($(cat $uuid_list))
#names=($(cat $net_names))

#uuid=${uuids[$SLURM_ARRAY_TASK_ID]}
#name=${names[$SLURM_ARRAY_TASK_ID]}

srun -l python $execdir/get_heat_matrix.py $OUTDIR $uuid $name $filter



