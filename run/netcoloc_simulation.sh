#!/bin/bash -l

DATADIR=/cellar/users/snwright/Data/NetColocTest/inputs/GO
execdir=/cellar/users/snwright/Git/rare_common/carva
OUTDIR=/cellar/users/snwright/Data/NetColocTest/outputs/overlap_test
mkdir -p $OUTDIR
uuid='d73d6357-e87b-11ee-9621-005056ae23aa'
name=pcnet2_0

traits=($(cat $DATADIR/trait_list_go_overlap_only_July11.txt))

#t=${traits[$SLURM_ARRAY_TASK_ID]}
t='GO:0000082_overlap40_relevance1.0_totalgenes150_repeat1_background0'
mkdir -p $OUTDIR

python $execdir/do_netcoloc.py --outdir $OUTDIR --indir $DATADIR --trait_rare $t --trait_common $t --uuid $uuid --net_name $name
