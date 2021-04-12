#!/usr/bin/env nextflow


////////////////////////////////////////////////////
/* --            Collect  parameters           -- */
////////////////////////////////////////////////////

Channel
    .fromPath("${params.proprdir}/*/propr_results.rds", checkIfExists:true)
    .map { it -> [ it.getParent().baseName, it ] }
    .combine( Channel.fromPath(params.keggfile) )
    .set{ ch_tokegg }


process kegg {

    tag "${subname}"
    publishDir "${params.outdir}/${subname}/kegg", mode: params.publish_dir_mode

    input:
    set val(subname), file(propr), file(kegg) from ch_tokegg

    output:
    set file('curve.jpg'), \
        file('pr.txt'), \
        file('roc.txt')

    script:
    """
    Rscript ${baseDir}/bin/kegg/kegg.R \
        --pro ${propr} \
        --kegg ${kegg}
    """
}

