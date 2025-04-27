#!/bin/bash -l
#SBATCH --job-name=netcolocX
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/Xnetcoloc_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/Xnetcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=00:30:00
#SBATCH --array=1-197%10

pair_file=$1
overlap_control=$2
DATADIR=$3
OUTDIR=$4
#DATADIR=/cellar/users/snwright/Data/RareCommon/inputs
execdir=/cellar/users/snwright/Git/rare_common/carva
#OUTDIR=/cellar/users/snwright/Data/RareCommon/outputs/netcoloc
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
name=pcnet2_0

readarray -t all_pairs < $DATADIR/$pair_file

pair=${all_pairs[$SLURM_ARRAY_TASK_ID]}

mkdir -p $OUTDIR
t=$(echo $pair | awk '{print $1}')
t2=$(echo $pair | awk '{print $2}')

srun -l python $execdir/do_netcoloc.py --outdir $OUTDIR --indir $DATADIR \
	--trait_rare $t --trait_common $t2 --uuid $uuid \
	--net_name $name --overlap_control $overlap_control

srun -l python $execdir/do_netcoloc.py --outdir $OUTDIR --indir $DATADIR \
	--trait_rare $t2 --trait_common $t --uuid $uuid \
	--net_name $name --overlap_control $overlap_control


