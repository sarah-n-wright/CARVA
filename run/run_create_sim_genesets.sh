#!/bin/bash -l
#SBATCH --job-name=genesets
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/genesets_%A.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/genesets_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500MB
#SBATCH --time=01:00:00

setfile=/cellar/users/snwright/Data/NetColocTest/Reference/go.test
total_genes=150
repeats=5
outdir=/cellar/users/snwright/Data/NetColocTest/inputs/GO/
nodefile=/cellar/users/snwright/Data/NetColocTest/Reference/pcnet2_0_nodelist.txt
execdir=/cellar/users/snwright/Git/rare_common/carva
# with an overlap of 50 there will be each gene set will have 100 genes, with 50 of them being mat
overlaps=( 0 30 60 ) # number of matching genes
relevance=( 1.0 0.6 0.3 )
#relevance=( 0 0.1 0.25 0.5 0.75 )
#relevance=( 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0 ) # persent of remaining genes to be from same gene set


for o in "${overlaps[@]}"
do
	for r in "${relevance[@]}"
	do
		# create the input files
		srun -l python $execdir/create_sim_genesets.py --setfile $setfile \
			--outdir $outdir --netnodefile $nodefile \
			--overlap $o --relevance $r \
			--totalgenes $total_genes --nrepeats $repeats --background degree
		## TODO
		#srun -l python do_netcoloc.py --outdir $outdir --indir $indir --trait_rare $suffix.1 \
		#	--trait_common $suffix.2 --uuid $uuid --net_name $net_name
	done
done
