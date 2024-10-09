#!/bin/bash -l
#SBATCH --job-name=netcolocX
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/Xnetcoloc_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/Xnetcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=00:30:00
#SBATCH --array=0-197%10

pair_file=$1
DATADIR=$2
OUTDIR=$3
#DATADIR=/cellar/users/snwright/Data/RareCommon/inputs
execdir=/cellar/users/snwright/Git/rare_common/carva
#OUTDIR=/cellar/users/snwright/Data/RareCommon/outputs/netcoloc
test_name=defaults
background=20000

readarray -t all_pairs < $DATADIR/$pair_file

pair=${all_pairs[$SLURM_ARRAY_TASK_ID]}

mkdir -p $OUTDIR
t=$(echo $pair | awk '{print $1}')
t2=$(echo $pair | awk '{print $2}')
echo $t
echo $t2

python $execdir/gene_overlap.py --datadir $DATADIR --raretrait $t --commontrait $t2 \
	--rare_th 1 --common_th 1 --min_genes 3 --outdir $OUTDIR \
	--test_name $test_name --background_N $background

python $execdir/gene_overlap.py --datadir $DATADIR --raretrait $t2 --commontrait $t \
	--rare_th 1 --common_th 1 --min_genes 3 --outdir $OUTDIR \
	--test_name $test_name --background_N $background

