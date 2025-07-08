#!/bin/bash -l
#SBATCH --job-name=paths
#SBATCH --output paths_%A_%a.out
#SBATCH --error paths_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48GB
#SBATCH --time=12:00:00
#SBATCH --array=0-3%2

uuids=('d73d6357-e87b-11ee-9621-005056ae23aa' '8b4b54fa-e87d-11ee-9621-005056ae23aa' '5f5da339-f14a-11ee-9621-005056ae23aa' 'db921c7d-f144-11ee-9621-005056ae23aa')
PWD=$(pwd)

outdir=$PWD/../outputs
outpref=('pcnet2' 'pcnet22' 'string' 'humannet')

net=${uuids[$SLURM_ARRAY_TASK_ID]}

pref=${outpref[$SLURM_ARRAY_TASK_ID]}

echo $pref $net

execdir=$PWD/../carva
srun -l python $execdir/find_shortest_paths.py --outdir $outdir --uuid $net --outpref $pref


