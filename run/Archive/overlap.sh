#!/bin/bash -l
#SBATCH --job-name=overlap
#SBATCH --output /cellar/users/snwright/Data/RareCommon/slurm/overlap_%A.out
#SBATCH --error /cellar/users/snwright/Data/RareCommon/slurm/overlap_%A.err
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00

DATADIR=/cellar/users/snwright/Data/NetColocTest/inputs/GO
execdir=/cellar/users/snwright/Git/rare_common/carva
OUTDIR=/cellar/users/snwright/Data/NetColocTest/outputs
mode=between # within or between
rare_th=1 #0.00001
common_th=1 #0.00000005
min_genes=3
test_name=defaults
background=20000 #9415

traits=($(cat $DATADIR/trait_pairs_cross_GO.txt))

parallel_path=~/anaconda3/bin/parallel

echo $test_name $mode $rare_th $common_th $min_genes >> $OUTDIR/overlap_params.txt
mkdir -p $OUTDIR/$test_name

process_overlap() {
    trait1=$1
    echo "TRAIT:" $trait1
    if [ $mode == "within" ]; then
        python $execdir/gene_overlap.py --datadir $DATADIR --raretrait $trait1 --commontrait $trait1 --outdir $OUTDIR/$test_name \
            --rare_th $rare_th --common_th $common_th --min_genes $min_genes --test_name $test_name --background_N $background
    else
        for trait2 in ${traits[@]}; do
            python $execdir/gene_overlap.py --datadir $DATADIR --raretrait $trait1 --commontrait $trait2 --outdir $OUTDIR/$test_name \
                --rare_th $rare_th --common_th $common_th --min_genes $min_genes --test_name $test_name --background_N $background
        done
    fi
}



export -f process_overlap
export DATADIR execdir OUTDIR mode traits rare_th common_th min_genes test_name background

$parallel_path -j 36 process_overlap ::: ${traits[@]}

cat $OUTDIR/$test_name/* > $OUTDIR/overlap_results.$test_name.txt

rm $OUTDIR/$test_name/*
rmdir $OUTDIR/$test_name
