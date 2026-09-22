#!/usr/bin/env nextflow

// Next steps: either get bioplex into Entrez, or fix NetColoc to allow specifying the column names

// Include modules
include { getNetColocMatrices } from './modules/getNetColocMatrices.nf'
include { geneOverlap } from './modules/geneOverlap.nf'
include { doNetColoc } from './modules/doNetColoc.nf'
include { combineColocFiles } from './modules/combineColocFiles.nf'
/*
* CARVAParameters
*/

params {
    network_uuid: String
    network_name: String
    network_filter: String = ''
    netprop_alpha: Float = 0.5
    OUTDIR: String
    common_traits_file: Path
    rare_traits_file: Path
    geneset_datadir: Path
    degree_match_binsize: Integer = 20
    gene_score_transform: String
    gene_score_normalization: String
    geneset_min_genes: Integer = 3
    netcoloc_overlap_control: String = 'bin'
    rare_gene_file_suff: String = '_RV'
    common_gene_file_suff: String = '_CV'
    netcoloc_individual_z_threshold: Float = 1.0
    netcoloc_combined_z_threshold: Float = 3.0
    output_suffix: String = ''
    stat_suffix: String = ''
    seed_gene_col: String = ''
    seed_score_col: String = ''
    batch_name: String = ''
}

/*
* Pipeline
*/

workflow {

    main:
    common_traits = channel.fromPath(params.common_traits_file)
                            .splitCsv(header:false)
                            .map { row -> row[0] }
    rare_traits = channel.fromPath(params.rare_traits_file)
                            .splitCsv(header:false)
                            .map { row -> row[0] }
    // Get NetColoc Matrices
    getNetColocMatrices(params.network_uuid, 
                        '.', 
                        params.network_name,
                        params.network_filter,
                        params.netprop_alpha as float)

    // Create gene sets
    geneOverlap(rare_traits, 
                common_traits, 
                '.', 
                params.geneset_datadir, 
                1, 1, 
                params.geneset_min_genes, 
                'gene_overlap', 
                20000)
    // Run network colocalization for each gene set and each network
    doNetColoc(rare_traits,
                common_traits,
                '.',
                params.geneset_datadir,
                getNetColocMatrices.out.netdir,
                params.degree_match_binsize, 
                params.network_uuid, 
                params.network_name, 
                params.gene_score_transform, 
                params.gene_score_normalization, 
                params.geneset_min_genes, 
                params.netcoloc_overlap_control, 
                params.rare_gene_file_suff, 
                params.common_gene_file_suff,
                params.output_suffix, 
                params.netcoloc_individual_z_threshold,
                params.netcoloc_combined_z_threshold,
                params.stat_suffix,
                params.seed_gene_col,
                params.seed_score_col)
    // Collate the results
    combineColocFiles(doNetColoc.out.coloc_results.collect(),
                        params.network_name,
                        params.batch_name)


    publish:
    // Network outputs
    w_prime = getNetColocMatrices.out.w_prime
    individual_heats = getNetColocMatrices.out.individual_heats
    nodes = getNetColocMatrices.out.nodes
    degrees = getNetColocMatrices.out.degrees
    // Gene set outputs - process line by line?
        // take datafile, and paired gene set identifiers. Rare common pair becomes the channel?
    overlap_results = geneOverlap.out.overlap_results
    // NetColoc outputs - takes identifiers
    rare_z_scores = doNetColoc.out.rare_z_scores
    common_z_scores = doNetColoc.out.common_z_scores
    coloc_results = doNetColoc.out.coloc_results
    // Combined results
    combined_results = combineColocFiles.out.combined_results

}

output {
// Emit outputs from getNetColocMatrices process
w_prime {
    path { "${params.OUTDIR}/inputs/" }
}

individual_heats {
    path { "${params.OUTDIR}/inputs/" }
}

nodes {
    path { "${params.OUTDIR}/inputs/" }
}

degrees{
    path { "${params.OUTDIR}/inputs/" }
}
// Gene set outputs
overlap_results {
    path { "${params.OUTDIR}/outputs/overlap/" }
}
// NetColoc outputs
rare_z_scores {
    path { "${params.OUTDIR}/outputs/z_scores/" }
}
common_z_scores {
    path { "${params.OUTDIR}/outputs/z_scores/" }
}
coloc_results {
    path { "${params.OUTDIR}/outputs/coloc_results/" }
}
// Combined results
combined_results {
    path { "${params.OUTDIR}/outputs/coloc_results/" }
}

}