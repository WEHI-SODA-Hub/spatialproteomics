#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WEHI-SODA-Hub/sp_segment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/WEHI-SODA-Hub/sp_segment
    Website: https://wehi-soda-hub.github.io/
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SP_SEGMENT  } from './workflows/sp_segment'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_sp_segment_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_sp_segment_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow WEHISODAHUB_SP_SEGMENT {

    take:
    samplesheet       // channel: samplesheet read in from --input
    aviti_samplesheet // channel: samplesheet read in from --aviti_input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    SP_SEGMENT (
        samplesheet,
        aviti_samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.aviti_input
    )

    //
    // WORKFLOW: Run main workflow
    //
    WEHISODAHUB_SP_SEGMENT (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.aviti_samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,

    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
