#!/bin/bash -l
#SBATCH --job-name=geneset_annot
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/geneset_annot_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/geneset_annot_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=2-00:00:00
#SBATCH --array=0

DATADIR=/cellar/users/snwright/Data/RareCommon/inputs/March_2025
execdir=/cellar/users/snwright/Git/rare_common/carva
OUTDIR=/cellar/users/snwright/Data/RareCommon/outputs/Features
geneset_file=not_best_rare.traitlist
geneset_file2=not_best_common.traitlist
#uuid_list=('d73d6357-e87b-11ee-9621-005056ae23aa' '8b4b54fa-e87d-11ee-9621-005056ae23aa' '5f5da339-f14a-11ee-9621-005056ae23aa' 'db921c7d-f144-11ee-9621-005056ae23aa')
net_files=('pcnet2_0.edgelist.tsv' 'pcnet2_2.edgelist.tsv' 'string.edgelist.tsv' 'humannet.edgelist.tsv')
netdir=/cellar/users/snwright/Data/RareCommon/Reference
name_list=('pcnet2_0' 'pcnet2_2' 'string' 'humannet')

#uuid=${uuid_list[$SLURM_ARRAY_TASK_ID]}
netpath=$netdir/${net_files[$SLURM_ARRAY_TASK_ID]}
name=${name_list[$SLURM_ARRAY_TASK_ID]}

# Pre-calculate network stats
#/usr/bin/time -v python $execdir/get_network_stats.py --uuid $uuid --net_name $name \
#		--outdir $OUTDIR


/usr/bin/time -v python $execdir/network_annotation2.py --indir $DATADIR \
		--outdir $OUTDIR/network_stats --update \
		--netfile $netpath \
		--geneset_list_file $geneset_file --net_name $name \
		--geneset_list_file2 $geneset_file2 \
		--clustering_file $OUTDIR/${name}_clustering_coefficients.csv \
		--paths_file $OUTDIR/${name}path_lengths.csv
