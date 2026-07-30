//
// Subworkflow with functionality specific to the WEHI-SODA-Hub/sp_segment pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version            // boolean: Display version and exit
    validate_params    // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args  //   array: List of positional nextflow CLI args
    outdir             //  string: The output directory where the results will be saved
    _input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Fail before any compute on a model that cannot be resolved.
    //
    // Built-in cellpose models, as reported by `cellpose.models.MODEL_NAMES` in
    // the pinned container (cellpose 4.2.1.1). Anything else given to
    // cellpose_pretrained_model is treated as a path to a custom model.
    //
    // cellpose_pretrained_model takes either a built-in model name or a path to
    // a custom model, so only validate as a path when it is not a known name.
    // Cellpose does not error on a --pretrained-model path that does not exist:
    // it falls back to its built-in weights, so a typo produced a full,
    // plausible, silently-wrong run.
    //
    def CELLPOSE_BUILTIN_MODELS = ['cpsam_v2', 'cpsam', 'cpdino', 'cpdino-vitb']
    if (params.cellpose_pretrained_model
        && !CELLPOSE_BUILTIN_MODELS.contains(params.cellpose_pretrained_model)
        && !file(params.cellpose_pretrained_model).exists()) {
        error(
            "cellpose_pretrained_model is neither a built-in model name nor an existing path: " +
            "${params.cellpose_pretrained_model}\n" +
            "Built-in models: ${CELLPOSE_BUILTIN_MODELS.join(', ')}"
        )
    }
    if (params.cellpose_models_dir && !file(params.cellpose_models_dir).exists()) {
        error("cellpose_models_dir does not exist: ${params.cellpose_models_dir}")
    }

    //
    // A model cache under /opt breaks the container, confusingly.
    //
    // Nextflow bind-mounts the host directory holding a staged input, so
    // --cellpose_models_dir /opt/... mounts the host's /opt over the
    // container's. The Cellpose image installs into /opt/conda, so that
    // shadows the interpreter and every task dies with
    // "sopa: command not found" -- an error that says nothing about the real
    // cause. /opt is a natural place to put a shared cache, so catch it here.
    //
    if (params.cellpose_models_dir && file(params.cellpose_models_dir).toAbsolutePath().toString().startsWith('/opt')) {
        error(
            "cellpose_models_dir must not be under /opt: ${params.cellpose_models_dir}\n" +
            "Nextflow bind-mounts the host directory containing a staged input, which would " +
            "mount the host's /opt over the container's and hide /opt/conda, where the Cellpose " +
            "image is installed. Tasks would fail with \"sopa: command not found\".\n" +
            "Put the cache somewhere else, e.g. /shared/cellpose_models."
        )
    }

    //
    // cellpose_model_type is retired rather than quietly ignored.
    //
    // Cellpose >=4.0.1 logs "model_type argument is not used in v4.0.1+" and
    // discards it, and sopa only reads it on the cellpose<4 path, so the old
    // 'cyto3' default silently produced cpsam results on every run. Failing is
    // better than repeating that.
    //
    if (params.cellpose_model_type) {
        error(
            "cellpose_model_type is no longer supported: cellpose 4 ignores --model-type.\n" +
            "Use --cellpose_pretrained_model instead " +
            "(${CELLPOSE_BUILTIN_MODELS.join(', ')}, or a path to a custom model)."
        )
    }

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications


    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {

    return input
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",


            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [


        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

