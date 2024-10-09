#!/bin/bash -l
#SBATCH --job-name=overlap
#SBATCH --output /cellar/users/snwright/Data/RareCommon/overlap_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/overlap_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00


DATADIR=/cellar/users/snwright/Data/RareCommon/

traits=(
    30630
    30640
    30830
    30840
    30880
    23110
    30050
    30140
    30250
    20015
    20151
    23098
    2453
    3063
    30680
    5133
    C3_PRIMARY_LYMPHOID_HEMATOPOIETIC
    C92
    D46
    D47
    XVII_MALFORMAT_ABNORMAL
)

files=(
    biomarkers-30630-both_sexes-irnt.tsv
    biomarkers-30640-both_sexes-irnt.tsv
    biomarkers-30830-both_sexes-irnt.tsv
    biomarkers-30840-both_sexes-irnt.tsv
    biomarkers-30880-both_sexes-irnt.tsv
    continuous-23110-both_sexes-irnt.tsv
    continuous-30050-both_sexes-irnt.tsv
    continuous-30140-both_sexes-irnt.tsv
    continuous-30250-both_sexes-irnt.tsv
    '20015_irnt.gwas.imputed_v3.both_sexes.tsv'
    '20151_irnt.gwas.imputed_v3.both_sexes.tsv'
    '23098_irnt.gwas.imputed_v3.both_sexes.tsv'
    '2453.gwas.imputed_v3.both_sexes.tsv'
    '3063_irnt.gwas.imputed_v3.both_sexes.tsv'
    '30680_irnt.gwas.imputed_v3.both_sexes.varorder.tsv'
    '5133_irnt.gwas.imputed_v3.both_sexes.tsv'
    C3_PRIMARY_LYMPHOID_HEMATOPOIETIC.gwas.imputed_v3.both_sexes.tsv
    C92.gwas.imputed_v3.both_sexes.tsv
    D46.gwas.imputed_v3.both_sexes.tsv
    D47.gwas.imputed_v3.both_sexes.tsv
    XVII_MALFORMAT_ABNORMAL.gwas.imputed_v3.both_sexes.tsv
)

pheno_type=(
    quant
    quant
    quant
    quant
    quant
    quant
    quant
    quant
    quant
    quant
    quant
    quant
    binary
    quant
    quant
    quant
    binary
    binary
    binary
    binary
    binary
)

#f=${files[$SLURM_ARRAY_TASK_ID]}
#t=${traits[$SLURM_ARRAY_TASK_ID]}
outFile=$DATADIR/pilot_overlap_results_matrix4.txt
echo -e "TraitC\tTraitR\tCommon\tRare\tOverlap\tOverlap_p" > $outFile

for c_id in {0..20}; do
    for r_id in {0..20}; do
        t1=${traits[$c_id]}
        t1_file=${files[$c_id]}
        t2=${traits[$r_id]}
        rare_table=${pheno_type[$r_id]}
        if [[ "$rare_table" == "binary" ]]; then
            rareData=$DATADIR/AZdata/Table16_binary_map.txt
            binary=1
        else
            rareData=$DATADIR/AZdata/Table17_quant_map.txt
            binary=0
        fi
        commonData=$DATADIR/NealeLab/pascal_outputs/$t1_file.sum.genescores.txt
        python /cellar/users/snwright/Git/rare_common/direct_overlap.py $commonData $rareData $t1 $t2 $binary $outFile
    done
done
