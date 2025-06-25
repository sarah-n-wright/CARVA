#!/bin/bash -l
#SBATCH --job-name=qnetcoloc
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/qnetcoloc_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/qnetcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

#DATADIR=/cellar/users/snwright/Data/RareCommon/inputs/testing
execdir=/cellar/users/snwright/Git/rare_common/carva
#OUTDIR=/cellar/users/snwright/Data/RareCommon/outputs/netcoloc/testing
#uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
#name=pcnet2_0
netdir=/cellar/users/snwright/Data/RareCommon/inputs/

config=$1
source run_configs/$1

traits1=($(cat $trait_list1))
traits2=($(cat $trait_list2))

t1=${traits1[$SLURM_ARRAY_TASK_ID]}
t2=${traits2[$SLURM_ARRAY_TASK_ID]}

file_list=${outdir}/${SLURM_ARRAY_JOB_ID}.files

echo $t1 $t2 $q $transform $normalization ${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID >> $trait_list.slurmIDs

mkdir -p $outdir

#echo -e "Trait_Common\tTrait_Rare\tNetwork\tMean_NPS\tNull_NPS\tp_NPS\tSize\tNull_Size\tp_Size" > $OUTDIR/netcoloc/pilot_netcoloc_results_$t_$t.txt

if [[ "$overlap" == 'bin' ]]; then
	overlap_control=bin
else
	overlap_control=remove
fi

if [[ $q -eq 1 ]]; then

/usr/bin/time -v srun -l python $execdir/do_carva_netcoloc.py --outdir $outdir \
	--indir $datadir --trait_rare $t1 --trait_common $t2 --netdir $netdir \
	--uuid $uuid --net_name $name --transform $transform --binsize 20 \
	--normalization $normalization --quant --min-genes 3 --overlap_control $overlap_control
echo qnetcoloc_${t1}_${t2}__q_${transform}_${normalization}.txt >> $file_list
/usr/bin/time -v srun -l python $execdir/do_carva_netcoloc.py --outdir $outdir \
	--indir $datadir --trait_rare $t2 --trait_common $t1 --netdir $netdir \
	--uuid $uuid --net_name $name --transform $transform --binsize 20 \
	--normalization $normalization --quant --min-genes 3 --overlap_control $overlap_control
echo qnetcoloc_${t2}_${t1}__q_${transform}_${normalization}.txt >> $file_list

else
/usr/bin/time -v srun -l python $execdir/do_carva_netcoloc.py --outdir $outdir \
	--indir $datadir --trait_rare $t1 --trait_common $t2 \
	--uuid $uuid --net_name $name --min-genes 3 --netdir $netdir
echo netcoloc_${t2}_${t1}__${normalization}.txt >> $file_list

/usr/bin/time -v srun -l python $execdir/do_carva_netcoloc.py --outdir $outdir \
	--indir $datadir --trait_rare $t2 --trait_common $t1 \
	--uuid $uuid --net_name $name --min-genes 3 --netdir $netdir
echo netcoloc_${t1}_${t2}__${normalization}.txt >> $file_list

fi
