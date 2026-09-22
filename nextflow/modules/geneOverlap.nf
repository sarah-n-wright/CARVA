process geneOverlap {

    conda 'netcoloc.yml'

    input:
        val rare_trait
        val common_trait
        val outdir
        val datadir
        val rare_th
        val common_th
        val min_genes
        val test_name
        val background

    output:
        path "${rare_trait}_${common_trait}.R_C_overlap.txt", emit: overlap_results

    script:
    """
    python ${baseDir}/../carva/gene_overlap.py --datadir ${datadir} \
    --raretrait ${rare_trait} --commontrait ${common_trait} \
	--rare_th ${rare_th} --common_th ${common_th} --min_genes ${min_genes} \
    --outdir ${outdir} \
	--test_name ${test_name} --background_N ${background}
    """

}