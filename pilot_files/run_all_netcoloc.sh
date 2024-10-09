#!/bin/bash -l
#SBATCH --job-name=netcoloc
#SBATCH --output /cellar/users/snwright/Data/RareCommon/netcoloc_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/netcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=03:00:00
#SBATCH --array=0-20%5

DATADIR="/cellar/users/snwright/Data/RareCommon"

traits=(
    30630
    30640
    30830
    30840
    30880
    23110
    30050
    30140
    30250
    20015
    20151
    23098
    2453
    3063
    30680
    5133
    C3_PRIMARY_LYMPHOID_HEMATOPOIETIC
    C92
    D46
    D47
    XVII_MALFORMAT_ABNORMAL
)

t=${traits[$SLURM_ARRAY_TASK_ID]}
commonSeeds=$DATADIR/common_seeds_$t.txt

echo -e "Trait_Common\tTrait_Rare\tMean_NPS\tNull_NPS\tp_NPS\tSize\tNull_Size\tp_Size" > $DATADIR/pilot_netcoloc_results_$t.txt

for t2 in ${traits[@]}; do
    rareSeeds=$DATADIR/rare_seeds_$t2.txt
    srun -l python do_netcoloc.py $DATADIR $commonSeeds $rareSeeds $t2 $t
done