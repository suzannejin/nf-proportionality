#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-proportionality
========================================================================================
 nf-proportionality Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/suzannejin/nf-proportionality.git
----------------------------------------------------------------------------------------
*/

log.info Headers.nf_core(workflow, params.monochrome_logs)


////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run main.nf --input 'rse_gene.Rdata' -profile singularity"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}


////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}


////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

// expression data
ch_data = Channel.fromPath(params.input, checkIfExists:true)

// tissue & donor options
def l_tissues = new File(params.tissues).text.readLines()
def l_donors = new File(params.donors).text.readLines()
Channel
    .fromList(l_tissues)
    .map { it -> ["--tissue", "'${it}'", it.replaceAll(" ", "_")] }
    .set { ch_tissues }
Channel
    .fromList(l_donors)
    .map { it -> ["--donor", "'${it}'", it] }
    .set { ch_donors }
if (params.test_propr){
    ch_tissues
        .mix(ch_donors)
        .set {ch_sub}
}else{
    Channel
        .value( ["", "", "allsample"])
        .mix(ch_tissues)
        .mix(ch_donors)
        .set {ch_sub}
}

// methods
Channel
    .from(params.methods)
    .map { it -> [it.tokenize('_')[0], it.tokenize('_')[1], it.tokenize('_')[2]] }
    .set{ ch_methods }

// if only run test
if (params.test_propr){
    ch_test = Channel.from("--test")
}else{
    ch_test = Channel.from("")
}

// other parameters
Channel
    .value( [params.cutoff_interval, params.interval_min, params.interval_max, params.permutation] )
    .set { ch_params }

// merge channels -> input channel
ch_data
    .combine(ch_sub)
    .combine(ch_methods)
    .combine(ch_test)
    .set{ch_input}


////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
// TODO nf-core: Report custom parameters here
summary['Input']            = params.input
summary['Tissues']          = params.tissues
summary["Donors"]           = params.donors
summary["Methods"]          = params.methods
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
summary['Config Files'] = workflow.configFiles.join(', ')

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-proportionality-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/proportionality Workflow Summary'
    section_href: 'https://github.com/nf-core/proportionality'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }


/*
 * Parse software version numbers
 */
process get_software_versions {

    publishDir "${params.tracedir}", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf('.csv') > 0) filename   // to only save the .csv file
                      else null
        }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file 'software_versions.csv'

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    R --version > v_R.txt
    Rscript -e "cat(noquote(as.character(packageVersion('recount'))))" > v_recount.txt
    Rscript -e "cat(noquote(as.character(packageVersion('zCompositions'))))" > v_zCompositions.txt
    Rscript -e "cat(noquote(as.character(packageVersion('propr'))))" > v_propr.txt
    Rscript -e "cat(noquote(as.character(packageVersion('edgeR'))))" > v_edgeR.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
 * STEP 1 - Preprocess count data: select, normalize and log transform
 */
process process_data {

    memory = { 4.GB * task.attempt }
    time = { 1.h * task.attempt }

    tag "${subname}"
    publishDir "${params.outdir}/${subname}", mode: params.publish_dir_mode,
        saveAs : { filename -> 
                        if ( filename == 'processed_data.Rdata' ) filename
                        else null
        }

    input:
    set file(input), \
        val(flag), \
        val(sub), \
        val(subname), \
        val(method), \
        val(norm), \
        val(ivar), \
        val(test) from ch_input

    output:
    set file('processed_data.rds'), val(method), val(subname) into ch_topropr

    script:
    """
    Rscript ${baseDir}/bin/propr/process_data.R \
        --data ${input} \
        --output processed_data.rds \
        ${flag} ${sub} \
        --norm ${norm} \
        --ivar ${ivar} \
        ${test}
    """
}

/*
 * STEP 2 - Run propr: compute association coefficients
 */
process propr {

    memory = { if ( subname == "allsample" ) 64.GB 
               else 12.GB * task.attempt
    }
    tag "${subname}"
    publishDir "${params.outdir}/${subname}", mode: params.publish_dir_mode

    input:
    set file(input), val(method), val(subname) from ch_topropr
    val(fdrparams) from ch_params

    output:
    file 'propr_results*'
    set val(subname), file('propr_results.rds') into ch_proprout

    script:
    """
    Rscript ${baseDir}/bin/propr/compute_propr.R \
        --data ${input} \
        --method ${method} \
        --cutoff_interval ${fdrparams[0]} \
        --interval_min ${fdrparams[1]} \
        --interval_max ${fdrparams[2]} \
        --permutation ${fdrparams[3]} 
    """
}


// create input for kegg process
ch_proprout
    .combine( Channel.fromPath(params.keggfile) )
    .combine( Channel.fromList(['-clique', '']))
    .set{ch_tokegg}

/*
 * STEP 3 - Compare coexpression matrix with KEGG pathways
 */
process kegg {

    memory = { 6.GB * task.attempt }
    time = { 30.min * task.attempt }

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



workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}
