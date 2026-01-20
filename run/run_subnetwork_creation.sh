#!/bin/bash -l
# Extract trait-specific subnetworks from colocalization results
# Creates subnetworks containing input genes and colocalized genes identified by NetColoc
# Usage: sbatch run_subnetwork_creation.sh <trait_list_file>
# Example: sbatch run_subnetwork_creation.sh significant_traits.tsv

#SBATCH --job-name=subnet
#SBATCH --output=subnet_%A.out
#SBATCH --error=subnet_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=03:00:00

# Input: TSV file with TraitR and TraitC columns listing trait pairs
traitlist_file=$1

# Set up directory paths
PWD=$(pwd)
DATADIR=$PWD/../outputs          # Input gene lists (<trait>_RV.txt and <trait>_CV.txt)
z_dir=$PWD/../outputs             # Z-score files from NetColoc analysis
execdir=$PWD/../carva
OUTDIR=$PWD/../outputs/subnetworks
uuid='d73d6357-e87b-11ee-9621-005056ae23aa' # PCNet 2.0

mkdir -p $OUTDIR

# Extract subnetworks using colocalization thresholds (zth=1, zzth=3)
python $execdir/create_subnetworks.py --network_uuid $uuid --trait_list_file $traitlist_file \
        --z_dir $z_dir --genelist_dir $DATADIR \
        --outputdir $OUTDIR --zth 1 --zzth 3 

