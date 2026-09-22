process combineColocFiles {

    input:
    path coloc_file
    val net_name
    val results_name

    output:
    path "${results_name}_${net_name}_coloc_results.tsv", emit: combined_results

    script:
    def header = "trait_rare\ttrait_common\tnetwork\ttransformation\tnormalization\tobserved_mean_NPS\tnull_mean_NPS\tp_mean_NPS\tobserved_size\tnull_size\tp_size"
    """
    echo "${header}" > ${results_name}_${net_name}_coloc_results.tsv
    cat ${coloc_file} >> ${results_name}_${net_name}_coloc_results.tsv
    """

}