#!/bin/bash -l
# Generate network propagation matrices for multiple biological networks
# Precomputes heat diffusion matrices required for NetColoc analysis
# Usage: sbatch get_netcoloc_matrices.sh

#SBATCH --job-name=matrix
#SBATCH --output matrix_%A_%a.out
#SBATCH --error matrix_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=01:00:00

PWD=$(pwd)
OUTDIR=$PWD/../outputs
execdir=$PWD/../carva

# Generate matrices for PCNet 2.0
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
name=pcnet2_0

srun -l python $execdir/get_heat_matrix.py --outdir $OUTDIR --uuid $uuid --name $name

# Generate matrices for PCNet 2.1
uuid='8b4b54fa-e87d-11ee-9621-005056ae23aa'
name=pcnet2_2

srun -l python $execdir/get_heat_matrix.py --outdir $OUTDIR --uuid $uuid --name $name

# Generate matrices for HumanNet
uuid='db921c7d-f144-11ee-9621-005056ae23aa'
name=humannet

srun -l python $execdir/get_heat_matrix.py --outdir $OUTDIR --uuid $uuid --name $name

# Generate matrices for STRING High Confidence (with edge filtering)
uuid='5f5da339-f14a-11ee-9621-005056ae23aa'
name=stringHC
filter=Score_700

srun -l python $execdir/get_heat_matrix.py --outdir $OUTDIR --uuid $uuid --name $name --filter $filter

# Generate matrices for OmniPath
uuid='1ed6be26-6bf3-11f0-a218-005056ae3c32'
name='omnipath'

srun -l python $execdir/get_heat_matrix.py --outdir $OUTDIR --uuid $uuid --name $name 




