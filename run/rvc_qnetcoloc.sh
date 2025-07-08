#!/bin/bash -l
#SBATCH --job-name=qnetcoloc
#SBATCH --output=qnetcoloc_%A_%a.out
#SBATCH --error=qnetcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=30G

PWD=$(pwd)
execdir=$PWD/../carva
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'  # PCNet 2.0
name=pcnet2_0
netdir=$PWD/../outputs/
raresuff='_RV'
commonsuff='_CV'

config=$1
source run_configs/$1

traitsR=($(cat $trait_listR))
traitsC=($(cat $trait_listC))

tR=${traitsR[$SLURM_ARRAY_TASK_ID]}
tC=${traitsC[$SLURM_ARRAY_TASK_ID]}

echo $tR $tC $q $transform $normalization ${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID >> $trait_list.slurmIDs

mkdir -p $outdir

file_list=${outdir}/${SLURM_ARRAY_JOB_ID}.files

if [[ "$overlap" == 'bin' ]]; then
	overlap_control=bin
else
	overlap_control=remove
fi


/usr/bin/time -v srun -l python $execdir/do_carva_netcoloc.py --outdir $outdir \
	--indir $datadir --trait_rare $tR --trait_common $tC \
	--netdir $netdir --binsize 20 \
	--uuid $uuid --net_name $name --transform $transform \
	--normalization $normalization --quant --min-genes 3 \
	--overlap_control $overlap_control --raresuff $raresuff --commonsuff $commonsuff

echo qnetcoloc_${tR}_${tC}__q_${transform}_${normalization}.txt >> $file_list
