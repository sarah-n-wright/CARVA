#!/bin/bash -l
#SBATCH --job-name=sumstats
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/sumstats_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/sumstats_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=01:00:00
#SBATCH --array=0-4,6-135%5

url_file=/cellar/users/snwright/Data/RareCommon/inputs/testing_gwas_sumstats_urls.txt

url_list=($(cat $url_file))

url=${url_list[$SLURM_ARRAY_TASK_ID]}

sh download_sumstats.sh $url
