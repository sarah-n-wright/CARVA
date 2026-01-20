#!/bin/bash -l
# Run Q-NetColoc analysis for rare and common variant gene sets
# Usage: sbatch rvc_qnetcoloc.sh <config_file>
# Example: sbatch rvc_qnetcoloc.sh config1.sh

#SBATCH --job-name=qnetcoloc
#SBATCH --output=qnetcoloc_%A_%a.out
#SBATCH --error=qnetcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --array=0-373%20

# Set up paths and network parameters
PWD=$(pwd)
execdir=$PWD/../carva
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'  # PCNet 2.0, modify for other networks.
name=pcnet2_0
netdir=$PWD/../outputs/
raresuff='_RV'
commonsuff='_CV'

# Load configuration file that defines: datadir, outdir, trait_listR, trait_listC, transform, normalization, overlap
config=$1 # This file must be formatted according to the template and be located in run/run_configs/
source run_configs/$1
traitsR=($(cat $trait_listR))
traitsC=($(cat $trait_listC))

# Select trait pair for this array task
tR=${traitsR[$SLURM_ARRAY_TASK_ID]}
tC=${traitsC[$SLURM_ARRAY_TASK_ID]}

# Log job information
echo $tR $tC $q $transform $normalization ${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID >> $trait_list.slurmIDs

mkdir -p $outdir

file_list=${outdir}/${SLURM_ARRAY_JOB_ID}.files

# Set overlap control method
if [[ "$overlap" == 'bin' ]]; then
	overlap_control=bin
else
	overlap_control=remove
fi

# Run Q-NetColoc with network propagation and permutation testing
# Note: Comment out the --quant flag for binary NetColoc
/usr/bin/time -v srun -l python $execdir/do_carva_netcoloc.py --outdir $outdir \
	--indir $datadir --trait_rare $tR --trait_common $tC \
	--netdir $netdir --binsize 20 \
	--uuid $uuid --net_name $name --transform $transform \
	--normalization $normalization \
	--min-genes 3 \
	--overlap_control $overlap_control --raresuff $raresuff --commonsuff $commonsuff \
	--quant

echo qnetcoloc_${tR}_${tC}__q_${transform}_${normalization}.txt >> $file_list
