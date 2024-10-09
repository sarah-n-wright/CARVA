#!/bin/bash -l
#SBATCH --job-name=netcoloc
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/netcoloc_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/netcoloc_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=500MB
#SBATCH --time=01:00:00
#SBATCH --array=0-1

DATADIR=/cellar/users/snwright/Data/RareCommon/inputs
execdir=/cellar/users/snwright/Git/rare_common
OUTDIR=/cellar/users/snwright/Data/RareCommon/outputs/netcoloc
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
name=pcnet2_0

traits=($(cat $DATADIR/overlap_traits_Jun27_min3_genes.txt))

t=${traits[$SLURM_ARRAY_TASK_ID]}

mkdir -p $OUTDIR

#echo -e "Trait_Common\tTrait_Rare\tNetwork\tMean_NPS\tNull_NPS\tp_NPS\tSize\tNull_Size\tp_Size" > $OUTDIR/netcoloc/pilot_netcoloc_results_$t_$t.txt

srun -l python $execdir/do_netcoloc.py --outdir $OUTDIR --indir $DATADIR --trait_rare $t --trait_common $t --uuid $uuid --net_name $name
