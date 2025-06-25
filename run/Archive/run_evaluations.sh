#!/bin/bash -l
#SBATCH --job-name=qnet_eval
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/q_eval_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/q_eval_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:10:00

#DATADIR=/cellar/users/snwright/Data/RareCommon/inputs/testing
execdir=/cellar/users/snwright/Git/rare_common/carva
#OUTDIR=/cellar/users/snwright/Data/RareCommon/outputs/netcoloc/testing

config=$1
source run_configs/$1

traits=($(cat $trait_list))

t=${traits[$SLURM_ARRAY_TASK_ID]}

echo $t $q $transform $normalization ${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID >> $trait_list.slurmIDs

if [[ $q -eq 1 ]]; then
	suffix=_q_${transform}_${normalization}
	outsuff=q_${normalization}
else
	suffix=_${transform}
	outsuff=binary
fi

/usr/bin/time -v srun -l python $execdir/ndcg_evaluation.py --genesetdir $datadir \
	--zscoredir $outdir --trait_study $t --outdir $outdir --zscoresuffix $suffix \
	--outsuffix $outsuff

