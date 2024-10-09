#!/bin/bash -l
#SBATCH --job-name=netcolocX
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/Xnetcoloc_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/Xnetcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500MB
#SBATCH --time=03:00:00
#SBATCH --array=0-100%10

DATADIR=/cellar/users/snwright/Data/RareCommon/inputs
execdir=/cellar/users/snwright/Git/rare_common/carva
OUTDIR=/cellar/users/snwright/Data/RareCommon/outputs/netcoloc
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
name=pcnet2_0

traits=($(cat $DATADIR/overlap_traits_Jun17_min3_genes.txt))

t=${traits[$SLURM_ARRAY_TASK_ID]}

mkdir -p $OUTDIR

#echo -e "Trait_Common\tTrait_Rare\tNetwork\tMean_NPS\tNull_NPS\tp_NPS\tSize\tNull_Size\tp_Size" > $OUTDIR/netcoloc/pilot_netcoloc_results_$t_$t.txt
for t2 in ${traits[@]}; do
	# check if t and t2 are the same
	if [ $t == $t2 ]; then
		continue
	else	
		srun -l python $execdir/do_netcoloc.py --outdir $OUTDIR --indir $DATADIR --trait_rare $t --trait_common $t2 --uuid $uuid --net_name $name
	fi
done
