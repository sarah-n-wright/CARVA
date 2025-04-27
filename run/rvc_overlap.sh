#!/bin/bash -l
#SBATCH --job-name=netcolocX
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/Xnetcoloc_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/Xnetcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=00:30:00
#SBATCH --array=0-372%10

rare_file=$1
common_file=$2
DATADIR=/cellar/users/snwright/Data/RareCommon/inputs/March_2025
OUTDIR=/cellar/users/snwright/Data/RareCommon/outputs/overlap/March_2025
#DATADIR=/cellar/users/snwright/Data/RareCommon/inputs
execdir=/cellar/users/snwright/Git/rare_common/carva
#OUTDIR=/cellar/users/snwright/Data/RareCommon/outputs/netcoloc
test_name=defaults
background=19000

rare_traits=($(cat $DATADIR/$rare_file))
common_traits=($(cat $DATADIR/$common_file))

tR=${rare_traits[$SLURM_ARRAY_TASK_ID]}
tC=${common_traits[$SLURM_ARRAY_TASK_ID]}

mkdir -p $OUTDIR
echo $tR
echo $tC

python $execdir/gene_overlap.py --datadir $DATADIR --raretrait $tR --commontrait $tC \
	--rare_th 1 --common_th 1 --min_genes 3 --outdir $OUTDIR \
	--test_name $test_name --background_N $background

