#!/bin/bash -l
# Calculate network topology statistics across multiple networks for rare and common variant gene sets
# Computes clustering, modularity, and assortativity for trait pairs across multiple networks
# Usage: sbatch annotations.sh

#SBATCH --job-name=geneset_annot
#SBATCH --output geneset_annot_%A_%a.out
#SBATCH --error geneset_annot_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=2-00:00:00
#SBATCH --array=0-3

# Set up paths and input files
PWD=$(pwd)
DATADIR=$PWD/../outputs
OUTDIR=$PWD/../outputs
execdir=$PWD/../carva

# Input: Two files listing trait names (one per line), matched by line number. These
# trait names should match the file prefixes of the common and rare variant gene set 
geneset_file=rare.traitlist       # File listing rare trait names (one per line)
geneset_file2=common.traitlist    # File listing common trait names (one per line)

# Network configurations (array task ID selects which network to use)
uuid_list=('d73d6357-e87b-11ee-9621-005056ae23aa' '8b4b54fa-e87d-11ee-9621-005056ae23aa' '5f5da339-f14a-11ee-9621-005056ae23aa' 'db921c7d-f144-11ee-9621-005056ae23aa')
name_list=('pcnet2_0' 'pcnet2_2' 'string' 'humannet')

uuid=${uuid_list[$SLURM_ARRAY_TASK_ID]}
name=${name_list[$SLURM_ARRAY_TASK_ID]}

# Pre-calculate clustering coefficients for all nodes (saved to prevent repeated calculation). Comment
# out if this has already been run. 
/usr/bin/time -v python $execdir/get_network_stats.py --uuid $uuid --net_name $name \
		--outdir $OUTDIR

# Calculate network statistics for all trait pairs
/usr/bin/time -v python $execdir/network_annotation.py --indir $DATADIR \
		--outdir $OUTDIR/network_stats --update \
		--uuid $uuid \
		--geneset_list_file $geneset_file --net_name $name \
		--geneset_list_file2 $geneset_file2 \
		--clustering_file $OUTDIR/${name}_clustering_coefficients.csv \
		--paths_file $OUTDIR/${name}path_lengths.csv
