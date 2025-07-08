#!/bin/bash -l
#SBATCH --job-name=matrix
#SBATCH --output matrix_%A_%a.out
#SBATCH --error matrix_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=01:00:00
#SBATCH --array=0

PWD=$(pwd)
OUTDIR=$PWD/../outputs
execdir=$PWD/../carva

## PCNET 2.0
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
name=pcnet2_0

srun -l python $execdir/get_heat_matrix.py $OUTDIR $uuid $name

## PCNET 2.1
uuid='8b4b54fa-e87d-11ee-9621-005056ae23aa'
name=pcnet2_2

srun -l python $execdir/get_heat_matrix.py $OUTDIR $uuid $name

## HumanNet
uuid='db921c7d-f144-11ee-9621-005056ae23aa'
name=humannet

srun -l python $execdir/get_heat_matrix.py $OUTDIR $uuid $name

## STRING High Confidence
uuid='5f5da339-f14a-11ee-9621-005056ae23aa'
name=stringHC
filter=Score_700

srun -l python $execdir/get_heat_matrix.py $OUTDIR $uuid $name $filter



