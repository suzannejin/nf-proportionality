#!/usr/bin/env nextflow


////////////////////////////////////////////////////
/* --            Collect  parameters           -- */
////////////////////////////////////////////////////
Channel
    .fromPath("${params.proprdir}/*/propr_results.rds", checkIfExists:true)
    .map { it -> [ it.getParent().baseName, it ] }
    .into { ch_proprout1; ch_proprout2 }

ch_proprout1
    .combine( Channel.fromPath(params.keggfile) )
    .combine( Channel.fromList(['clique', 'NA']) )  // kegg pathways as clique
    .combine( Channel.fromList([10, 2]) )           // min kegg pathway size
    .combine( Channel.fromList([1000, 100, 'NA'])  )     // max kegg pathway size
    .set{ ch_tokegg }

ch_proprout2
    .combine( Channel.fromList([10, 2]) )           // min GO size
    .combine( Channel.fromList([100, 'NA'])  )      // max GO size
    .set{ ch_togo }


process kegg {

    tag "${subname}"
    publishDir "${params.outdir}/${subname}/kegg/kegg-${clique}-minK${minK}-maxK${maxK}", mode: params.publish_dir_mode

    input:
    set val(subname), \
        file(propr), \
        file(kegg), \
        val(clique), \
        val(minK), \
        val(maxK) from ch_tokegg

    output:
    set file('curve.jpg'), \
        file('pr.txt'), \
        file('roc.txt')

    script:
    def clique_var = clique == 'clique' ? '--clique' : ''
    def maxK_var = maxK == 'NA' ? '' : "--maxK ${maxK}"

    """
    Rscript ${baseDir}/bin/kegg/kegg.R \
        --pro ${propr} \
        --kegg ${kegg} \
        ${clique_var} \
        --minK ${minK} \
        ${maxK_var}
    """
}

process go {

    tag "${subname}"
    publishDir "${params.outdir}/${subname}/go/go-minK${minK}-maxK${maxK}", mode: params.publish_dir_mode

    input:
    set val(subname), \
        file(propr), \
        val(minK), \
        val(maxK) from ch_togo

    output:
    set file('curve.jpg'), \
        file('pr.txt'), \
        file('roc.txt')

    script:
    def maxK_var = maxK == 'NA' ? '' : "--maxK ${maxK}"

    """
    Rscript ${baseDir}/bin/go/go.R \
        --pro ${propr} \
        --minK ${minK} \
        ${maxK_var}
    """
}
