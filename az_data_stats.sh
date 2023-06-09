DATADIR=/cellar/users/snwright/Data/RareCommon/AZdata

echo "#Phenotypes" > $DATADIR/num_phenotypes.txt
echo "#Genes" > $DATADIR/num_genes.txt
echo "#Unique_Combos" > $DATADIR/num_combos.txt
echo "Dataset" > $DATADIR/index.txt
echo "AtLeast2Genes" > $DATADIR/two_genes_plus.txt

## TABLE 12
f12=$DATADIR/Table12_onc.txt
echo $(wc -l $f12) >> $DATADIR/index.txt 
cut -f4 $f12 | sort | uniq | wc -l >> $DATADIR/num_phenotypes.txt 
cut -f5 $f12 | sort | uniq | wc -l >> $DATADIR/num_genes.txt 
cut -f4,5 $f12 | sort | uniq | wc -l >> $DATADIR/num_combos.txt 
cut -f4,5 $f12 | sort -k 1 | uniq | datamash -s -g 1 count 2  > $f12.pheno_counts
cat $f12.pheno_counts | datamash -s -g 2 count 1 | sort -n -k 1 > $f12.num_phenos 


## TABLE 8 - binary
f8b=$DATADIR/Table8_top_binary.txt
echo $(wc -l $f8b) >> $DATADIR/index.txt 
cut -f2 $f8b | sort | uniq | wc -l >> $DATADIR/num_phenotypes.txt 
cut -f4 $f8b | sort | uniq | wc -l >> $DATADIR/num_genes.txt 
cut -f2,4 $f8b | sort | uniq | wc -l >> $DATADIR/num_combos.txt 
cut -f2,4 $f8b | sort -k 1 | uniq | datamash -s -g 1 count 2  > $f8b.pheno_counts
cat $f8b.pheno_counts | datamash -s -g 2 count 1 | sort -n -k 1 > $f8b.num_phenos 

## TABLE 8 - quantitative
f8q=$DATADIR/Table8_top_quant.txt
echo $(wc -l $f8q) >> $DATADIR/index.txt 
cut -f5 $f8q | sort | uniq | wc -l >> $DATADIR/num_phenotypes.txt 
cut -f2 $f8q | sort | uniq | wc -l >> $DATADIR/num_genes.txt
cut -f2,5 $f8q | sort | uniq | wc -l >> $DATADIR/num_combos.txt
cut -f2,5 $f8q | sort -k 2 | uniq | datamash -s -g 2 count 1 > $f8q.pheno_counts
cat $f8q.pheno_counts | datamash -s -g 2 count 1 | sort -n -k 1 > $f8q.num_phenos 


## Table 17 - pan quant
f17=$DATADIR/Table17_pan_quant.txt
echo $(wc -l $f17) >> $DATADIR/index.txt 
cut -f3 $f17 | sort | uniq | wc -l >> $DATADIR/num_phenotypes.txt 
cut -f1 $f17 | sort | uniq | wc -l >> $DATADIR/num_genes.txt 
cut -f1,3 $f17 | sort | uniq | wc -l >> $DATADIR/num_combos.txt 
cut -f1,3 $f17 | sort -k 2 | uniq | datamash -s -g 2 count 1 > $f17.pheno_counts
cat $f17.pheno_counts | datamash -s -g 2 count 1 | sort -n -k 1 > $f17.num_phenos 

## Table 16 - binary pan
f16=$DATADIR/Table16_pan_binary.txt
echo $(wc -l $f16) >> $DATADIR/index.txt 
cut -f1 $f16 | sort | uniq | wc -l >> $DATADIR/num_phenotypes.txt 
cut -f2 $f16 | sort | uniq | wc -l >> $DATADIR/num_genes.txt 
cut -f1,2 $f16 | sort | uniq | wc -l >> $DATADIR/num_combos.txt 
cut -f1,2 $f16 | sort -k 1 | uniq | datamash -s -g 1 count 2 > $f16.pheno_counts
cat $f16.pheno_counts | datamash -s -g 2 count 1 | sort -n -k 1 > $f16.num_phenos 

paste $DATADIR/index.txt $DATADIR/num_phenotypes.txt $DATADIR/two_genes_plus.txt $DATADIR/num_genes.txt $DATADIR/num_combos.txt > $DATADIR/summary.txt

