#!/usr/bin/env nextflow


////////////////////////////////////////////////////
/* --            Collect  parameters           -- */
////////////////////////////////////////////////////

Channel
    .fromPath("${params.proprdir}/*/propr_results.rds", checkIfExists:true)
    .map { it -> [ it.getParent().baseName, it ] }
    .combine( Channel.fromPath(params.keggfile) )
    .combine( Channel.fromList(['-clique', '']))
    .set{ ch_tokegg }


process kegg {

    tag "${subname}"
    publishDir "${params.outdir}/${subname}/kegg${clique}", mode: params.publish_dir_mode

    input:
    set val(subname), file(propr), file(kegg), val(clique) from ch_tokegg

    output:
    set file('curve.jpg'), \
        file('pr.txt'), \
        file('roc.txt')

    script:
    def clique_var = clique == '-clique' ? '--clique' : ''

    """
    Rscript ${baseDir}/bin/kegg/kegg.R \
        --pro ${propr} \
        --kegg ${kegg} \
        ${clique_var}
    """
}

