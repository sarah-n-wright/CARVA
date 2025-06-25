#!/bin/bash -l
#SBATCH --job-name=qnetcoloc
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/qnetcoloc_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/qnetcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=32G

#DATADIR=/cellar/users/snwright/Data/RareCommon/inputs/testing
execdir=/cellar/users/snwright/Git/rare_common/carva
outdir=/cellar/users/snwright/Data/RareCommon/outputs/netcoloc/debug
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
name=pcnet2_0
netdir=/cellar/users/snwright/Data/RareCommon/inputs/



traitsR=($(cat $trait_listR))
traitsC=($(cat $trait_listC))

tR=36280733_EFO_0004574 #34375979_EFO_0004533
tC=GCST90025953_EFO_0004574 #GCST90025947_EFO_0004533

36280733_EFO_0004574_GCST90025953_EFO_0004574

mkdir -p $outdir

#echo -e "Trait_Common\tTrait_Rare\tNetwork\tMean_NPS\tNull_NPS\tp_NPS\tSize\tNull_Size\tp_Size" > $OUTDIR/netcoloc/pilot_netcoloc_results_$t_$t.txt


overlap_control=bin

/usr/bin/time -v srun -l python $execdir/do_carva_netcoloc.py --outdir $outdir \
	--indir $netdir/March_2025 --trait_rare $tR --trait_common $tC \
	--netdir $netdir --binsize 20 \
	--uuid $uuid --net_name $name --transform neglog10 \
	--normalization sum --quant --min-genes 3 \
	--overlap_control $overlap_control

