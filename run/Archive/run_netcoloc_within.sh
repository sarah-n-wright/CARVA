#!/bin/bash -l
#SBATCH --job-name=netcoloc
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/netcoloc_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/netcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=01:00:00
#SBATCH --array=11,13,18%2

DATADIR=/cellar/users/snwright/Data/RareCommon/inputs/testing
execdir=/cellar/users/snwright/Git/rare_common/carva
OUTDIR=/cellar/users/snwright/Data/RareCommon/outputs/netcoloc/testing
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
name=pcnet2_0

traits=($(cat $DATADIR/nearestGeneTH8.traitlist))

t=${traits[$SLURM_ARRAY_TASK_ID]}

mkdir -p $OUTDIR

#echo -e "Trait_Common\tTrait_Rare\tNetwork\tMean_NPS\tNull_NPS\tp_NPS\tSize\tNull_Size\tp_Size" > $OUTDIR/netcoloc/pilot_netcoloc_results_$t_$t.txt

/usr/bin/time -v srun -l python $execdir/do_netcoloc.py --outdir $OUTDIR --indir $DATADIR --trait_rare $t --trait_common $t --uuid $uuid --net_name $name

