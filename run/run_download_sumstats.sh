#!/bin/bash -l
#SBATCH --job-name=sumstats
#SBATCH --output sumstats_%A_%a.out
#SBATCH --error sumstats_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=01:00:00
#SBATCH --array=0-135%5

PWD=$(pwd)
OUTDIR=$PWD/../outputs
INDIR=$PWD/../Reference_Data

url_file=$INDIR/testing_gwas_sumstats_urls.txt

url_list=($(cat $url_file))

url=${url_list[$SLURM_ARRAY_TASK_ID]}


sh download_sumstats.sh $url $OUTDIR
