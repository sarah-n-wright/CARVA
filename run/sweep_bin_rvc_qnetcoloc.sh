#!/bin/bash -l
#SBATCH --job-name=qnetcoloc
#SBATCH --output=qnetcoloc_%A_%a.out
#SBATCH --error=qnetcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --array=0-50%5

PWD=$(pwd)
execdir=$PWD/../carva
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
name=pcnet2_0
netdir=$PWD/../outputs/

config=$1
source run_configs/$1

traitsR=($(cat $trait_list1))
traitsC=($(cat $trait_list2))

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

bin_sizes=(5 10 20 40 80 160)
for bin in ${bin_sizes[@]}; do
	test_suff=bin$bin

	/usr/bin/time -v srun -l python $execdir/do_carva_netcoloc.py --outdir $outdir \
		--indir $datadir --trait_rare $tR --trait_common $tC \
		--netdir $netdir --suffix $test_suff \
		--uuid $uuid --net_name $name --transform $transform \
		--normalization $normalization --quant --min-genes 3 \
		--overlap_control $overlap_control --binsize $bin

	echo qnetcoloc_${tR}_${tC}__q_${transform}_${normalization}.txt >> $file_list

done
