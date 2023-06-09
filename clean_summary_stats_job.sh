#!/bin/bash -l
#SBATCH --job-name=sum_stats
#SBATCH --output /cellar/users/snwright/Data/RareCommon/stats_%A.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/stats_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1:00:00
#SBATCH --array=0-5
file_only=$1
# Make sure rsid map files are unzipped before running

# files=(
#     biomarkers-30640-both_sexes-irnt.tsv
#     biomarkers-30830-both_sexes-irnt.tsv
#     biomarkers-30840-both_sexes-irnt.tsv
#     biomarkers-30880-both_sexes-irnt.tsv
#     continuous-23110-both_sexes-irnt.tsv
#     continuous-30050-both_sexes-irnt.tsv
#     continuous-30140-both_sexes-irnt.tsv
#     continuous-30250-both_sexes-irnt.tsv
#     icd10-E80-both_sexes.tsv
# )

# files=(
#     '2453.gwas.imputed_v3.both_sexes.tsv'
#     C3_PRIMARY_LYMPHOID_HEMATOPOIETIC.gwas.imputed_v3.both_sexes.tsv
#     C92.gwas.imputed_v3.both_sexes.tsv
#     D46.gwas.imputed_v3.both_sexes.tsv
#     D47.gwas.imputed_v3.both_sexes.tsv
#     XVII_MALFORMAT_ABNORMAL.gwas.imputed_v3.both_sexes.tsv
# )

files=(
    '20015_irnt.gwas.imputed_v3.both_sexes.tsv'
    '20151_irnt.gwas.imputed_v3.both_sexes.tsv'
    '23098_irnt.gwas.imputed_v3.both_sexes.tsv'
    '3063_irnt.gwas.imputed_v3.both_sexes.tsv'
    '5133_irnt.gwas.imputed_v3.both_sexes.tsv'
    '30680_irnt.gwas.imputed_v3.both_sexes.varorder.tsv'
)

f=${files[$SLURM_ARRAY_TASK_ID]}

sh /cellar/users/snwright/Git/rare_common/clean_summary_stats.sh $f $file_only
