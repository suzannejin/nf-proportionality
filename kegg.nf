#!/usr/bin/env nextflow


////////////////////////////////////////////////////
/* --            Collect  parameters           -- */
////////////////////////////////////////////////////
Channel
    .fromPath("${baseDir}/results_enzyme/*/*/propr_results.rds", checkIfExists:true)
    .map { it -> [ it.getParent().baseName, it.getParent(), it ] }
    .into { ch_proprout1; ch_proprout2; ch_proprout3 }

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

ch_proprout3
    .combine( Channel.fromPath(params.hippiefile) )
    .combine( Channel.fromList([0,1,2,3]) )           // HIPPIE score confidence level
    .set{ ch_tohippie }


process kegg {

    tag "${subname}"
    publishDir "${maindir}/kegg/kegg-${clique}-minK${minK}-maxK${maxK}", mode: params.publish_dir_mode

    input:
    set val(subname), \
        val(maindir), \
        file(propr), \
        file(kegg), \
        val(clique), \
        val(minK), \
        val(maxK) from ch_tokegg

    output:
    set file('curve.jpg'), \
        file('pr.txt'), \
        file('roc.txt')

    when:
        !file("${maindir}/kegg/kegg-${clique}-minK${minK}-maxK${maxK}/pr.txt").exists()

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
    publishDir "${maindir}/go/go-minK${minK}-maxK${maxK}", mode: params.publish_dir_mode

    input:
    set val(subname), \
        val(maindir), \
        file(propr), \
        val(minK), \
        val(maxK) from ch_togo

    output:
    set file('curve.jpg'), \
        file('pr.txt'), \
        file('roc.txt')

    when:
    !file("${maindir}/go/go-minK${minK}-maxK${maxK}/pr.txt").exists()

    script:
    def maxK_var = maxK == 'NA' ? '' : "--maxK ${maxK}"

    """
    Rscript ${baseDir}/bin/go/go.R \
        --pro ${propr} \
        --minK ${minK} \
        ${maxK_var}
    """
}


process hippie {

    tag "${subname}"
    publishDir "${maindir}/hippie/hippie-confidence${confidence}", mode: params.publish_dir_mode

    input:
    set val(subname), \
        val(maindir), \
        file(propr), \
        file(hippie), \
        val(confidence) from ch_tohippie

    output:
    set file('curve.jpg'), \
        file('pr.txt'), \
        file('roc.txt')

    when:
    !file("${maindir}/hippie/hippie-confidence${confidence}/pr.txt").exists()

    script:
    """
    Rscript ${baseDir}/bin/hippie/hippie.R \
        --pro ${propr} \
        --hippie ${hippie} \
        --confidence ${confidence}
    """
}

