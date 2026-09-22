process doNetColoc {

    //errorStrategy 'ignore'
    conda 'netcoloc.yml'

    input:
        val rare_trait
        val common_trait
        val outdir
        val datadir
        path netdir
        val binsize
        val uuid
        val net_name
        val transform
        val normalization
        val min_genes
        val overlap_control
        val rare_gene_file_suff
        val common_gene_file_suff
        val output_suffix
        val zcoloc
        val z1z2
        val stat_suffix
        val seed_gene_col
        val seed_score_col

    output:
        path "${rare_trait}_z_RV_q_${transform}_${normalization}.tsv", optional: true, emit: rare_z_scores
        path "${common_trait}_z_CV_q_${transform}_${normalization}.tsv", optional: true, emit: common_z_scores
        path "qnetcoloc_${rare_trait}_${common_trait}__q_${transform}_${normalization}.txt", optional: true, emit: coloc_results  

    script:
    def suffix_arg = output_suffix ? "--output_suffix ${output_suffix}" : ""
    def stat_arg = stat_suffix ? "--stat_suffix ${stat_suffix}" : ""
    def seed_gene_col_arg = seed_gene_col ? "--seed_gene_col \"${seed_gene_col}\"" : ""
    def seed_score_col_arg = seed_score_col ? "--seed_score_col \"${seed_score_col}\"" : ""

    """
    python ${baseDir}/../carva/do_carva_netcoloc.py --outdir ${outdir} \
	--indir ${datadir} --trait_rare ${rare_trait} --trait_common ${common_trait} \
	--netdir ${netdir} --binsize ${binsize} \
	--uuid ${uuid} --net_name ${net_name} --transform ${transform} \
	--normalization ${normalization} \
	--min-genes ${min_genes} --overlap_control ${overlap_control} --raresuff ${rare_gene_file_suff} \
    --commonsuff ${common_gene_file_suff} --quant --zcoloc ${zcoloc} --z1z2 ${z1z2} \
    ${suffix_arg} ${stat_arg} ${seed_gene_col_arg} ${seed_score_col_arg}
    """

}