/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATEREPORT           } from '../../modules/local/createreport/main'

workflow ANALYSE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    report_template_ch = "${moduleDir}/../../templates/report_template.qmd"
    CREATEREPORT(
        ch_samplesheet,
        report_template_ch
    )

    emit:
    rds      = CREATEREPORT.out.rds           // channel: [ val(meta), path(html), path(files) ]
    report   = CREATEREPORT.out.report        // channel: [ val(meta), path(rds) ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

