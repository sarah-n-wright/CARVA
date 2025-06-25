#!/bin/bash -l
#SBATCH --job-name=paths
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/paths_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/paths_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48GB
#SBATCH --time=12:00:00
#SBATCH --array=0-3%2

#uuids=('d73d6357-e87b-11ee-9621-005056ae23aa' '8b4b54fa-e87d-11ee-9621-005056ae23aa' '5f5da339-f14a-11ee-9621-005056ae23aa' 'db921c7d-f144-11ee-9621-005056ae23aa')
net_files=('pcnet2_0.edgelist.tsv' 'pcnet2_2.edgelist.tsv' 'string.edgelist.tsv' 'humannet.edgelist.tsv')
netdir=/cellar/users/snwright/Data/RareCommon/Reference

outdir=/cellar/users/snwright/Data/RareCommon/outputs/Features/
outpref=('pcnet2' 'pcnet22' 'string' 'humannet')

#net=${uuids[$SLURM_ARRAY_TASK_ID]}
net_input=$netdir/${net_files[$SLURM_ARRAY_TASK_ID]}

pref=${outpref[$SLURM_ARRAY_TASK_ID]}

echo $pref $net

execdir=/cellar/users/snwright/Git/rare_common/carva
#srun -l python $execdir/find_shortest_paths.py --outdir $outdir --uuid $net --outpref $pref
srun -l python $execdir/find_shortest_paths.py --outdir $outdir --netfile $net_input --outpref $pref


