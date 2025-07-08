#!/bin/bash -l
#SBATCH --job-name=subnet
#SBATCH --output=subnet_%A.out
#SBATCH --error=subnet_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=03:00:00

traitlist_file=$1
PWD=$(pwd)

DATADIR=$PWD/../outputs
z_dir=$PWD/../outputs
execdir=$PWD/../carva
OUTDIR=$PWD/../outputs/subnetworks
uuid='d73d6357-e87b-11ee-9621-005056ae23aa' # PCNet 2.0

echo 'TEST'
mkdir -p $OUTDIR

python $execdir/create_subnetworks.py --network_uuid $uuid --trait_list_file $traitlist_file \
        --z_dir $z_dir --genelist_dir $DATADIR \
        --outputdir $OUTDIR --zth 1 --zzth 3 

