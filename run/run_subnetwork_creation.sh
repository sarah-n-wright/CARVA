#!/bin/bash -l
#SBATCH --job-name=subnet
#SBATCH --output=/cellar/users/snwright/Data/RareCommon/slurm/subnet_%A.out
#SBATCH --error=/cellar/users/snwright/Data/RareCommon/slurm/subnet_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=03:00:00

traitlist_file=/cellar/users/snwright/Data/RareCommon/outputs/netcoloc/March_2025/coloc_repeat_best_results.tsv

DATADIR=/cellar/users/snwright/Data/RareCommon/inputs/March_2025
z_dir=/cellar/users/snwright/Data/RareCommon/outputs/netcoloc/March_2025
execdir=/cellar/users/snwright/Git/rare_common/carva
OUTDIR=/cellar/users/snwright/Data/RareCommon/outputs/subnetworks/pcnet2_0
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'

echo 'TEST'
mkdir -p $OUTDIR


python $execdir/create_subnetworks.py --network_uuid $uuid --trait_list_file $traitlist_file \
        --z_dir $z_dir --genelist_dir $DATADIR \
        --outputdir $OUTDIR --zth 1 --zzth 3 

