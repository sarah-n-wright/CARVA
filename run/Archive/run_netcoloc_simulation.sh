#!/bin/bash -l
#SBATCH --job-name=netcoloctest
#SBATCH --output /cellar/users/snwright/Data/NetColocTest/slurm/netcoloc_%A.out
#SBATCH --error /cellar/users/snwright/Data/NetColocTest/slurm/netcoloc_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=00:30:00
#SBATCH --array=0-1184%12

trait_file=$1
overlap_control=$2
out_folder=$3
DATADIR=/cellar/users/snwright/Data/NetColocTest/inputs/GO
execdir=/cellar/users/snwright/Git/rare_common/carva
OUTDIR=/cellar/users/snwright/Data/NetColocTest/outputs/$out_folder
mkdir -p $OUTDIR
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
name=pcnet2_0

traits=($(cat $DATADIR/$1))

t=${traits[$SLURM_ARRAY_TASK_ID]}

mkdir -p $OUTDIR

python $execdir/do_netcoloc.py --outdir $OUTDIR --indir $DATADIR --trait_rare $t \
	--trait_common $t --uuid $uuid --net_name $name --overlap_control $overlap_control
