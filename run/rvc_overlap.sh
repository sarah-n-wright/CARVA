#!/bin/bash -l
# Calculate gene overlap significance between paired lists of rare and common variant gene sets
# Usage: sbatch rvc_overlap.sh <rare_traits_file> <common_traits_file>
# Example: sbatch rvc_overlap.sh rare_trait_list.txt common_trait_list.txt

#SBATCH --job-name=netcolocX
#SBATCH --output Xnetcoloc_%A_%a.out
#SBATCH --error Xnetcoloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=00:30:00
#SBATCH --array=0-373 Run 373 tasks (adjust to match number of trait pairs)

# Input: Two files listing trait names (one per line), matched by line number. These
# trait names should match the file prefixes of the common and rare variant gene set 
# <trait>_RV.txt and <trait>_CV.txt files.
rare_file=$1
common_file=$2

# Set up directory paths
PWD=$(pwd)
DATADIR=$PWD/../outputs  # Input directory containing <trait>_RV.txt and <trait>_CV.txt files
OUTDIR=$PWD/../outputs   # Output directory for overlap statistics
execdir=$PWD/../carva    # Location of CARVA scripts

# Analysis parameters
test_name=defaults       # Suffix for output files
background=19000         # Background gene set size for hypergeometric test

# Load trait lists and select pair for this array task
rare_traits=($(cat $DATADIR/$rare_file))
common_traits=($(cat $DATADIR/$common_file))

tR=${rare_traits[$SLURM_ARRAY_TASK_ID]}
tC=${common_traits[$SLURM_ARRAY_TASK_ID]}

mkdir -p $OUTDIR
echo $tR
echo $tC

# Run hypergeometric test for gene overlap between rare and common variant sets
# Note rare_th=1 and common_th=1 results in no additonal filtering of the gene sets
python $execdir/gene_overlap.py --datadir $DATADIR --raretrait $tR --commontrait $tC \
	--rare_th 1 --common_th 1 --min_genes 3 --outdir $OUTDIR \
	--test_name $test_name --background_N $background

