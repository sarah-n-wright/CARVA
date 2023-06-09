##TODO - only do second join on missing snps. 

DATADIR=/cellar/users/snwright/Data/RareCommon/NealeLab
file=$1
just_produce_output=$2
rsid_map1=$DATADIR/full_rsid_map_sort_1.txt
rsid_map3=$DATADIR/full_rsid_map_sort_3.txt

if [ $just_produce_output -eq 0 ]
then
    #mv $DATADIR/$file.bgz $DATADIR/$file.gz
    #gunzip $DATADIR/$file.gz
    #echo "Unzipped $file"
# to update to use config file
    #cut -f1,2,3,4,8 $DATADIR/$file | awk -F "\t" '(NR>1){print $0 "\t" $1":"$2":"$3":"$4 }' | sort -k 6 > $DATADIR/$file.subset
    #echo "Created SNP ids for $file"
    sort -k 1 $DATADIR/$file > $DATADIR/$file.sort
# convert the snp ids
    join -1 1 -2 1 $rsid_map1 $DATADIR/$file.sort > $DATADIR/$file.snp1
    echo "Finished first rsid conversion for $file"
    join -1 3 -2 1 $rsid_map3 $DATADIR/$file.sort > $DATADIR/$file.snp3
    echo "Finished second rsid conversion for $file"

else
    #gzip $DATADIR/$file
    #TODO fix this because different files have different columns
    awk -v outfile=$DATADIR/$file.pascal_input '{print $2 "\t" $13 > outfile}' $DATADIR/$file.snp1
    awk -v outfile=$DATADIR/$file.pascal_input '{print $3 "\t" $13 >> outfile}' $DATADIR/$file.snp3

    echo "Created pascal input for $file"
fi

echo "Finished $file"
