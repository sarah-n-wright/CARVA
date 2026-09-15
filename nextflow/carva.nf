#!/usr/bin/env nextflow

// Include modules
include { getNetColocMatrices } from './modules/getNetColocMatrices.nf'
/*
* CARVAParameters
*/

params {
    network_uuid: String
    network_name: String
    network_filter: String = ''
    netprop_alpha: Float = 0.5
    OUTDIR: String
}

/*
* Pipeline
*/

workflow {

    main:
    // Get NetColoc Matrices
    getNetColocMatrices(params.network_uuid, 
                        '.', 
                        params.network_name,
                        params.network_filter,
                        params.netprop_alpha as float)

    publish:
    w_prime = getNetColocMatrices.out.w_prime
    individual_heats = getNetColocMatrices.out.individual_heats
    nodes = getNetColocMatrices.out.nodes
    degrees = getNetColocMatrices.out.degrees
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

}