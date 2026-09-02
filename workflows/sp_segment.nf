/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryMap          } from 'plugin/nf-schema'
include { BACKGROUNDSUBTRACT        } from '../subworkflows/local/backgroundsubtract'
include { MESMER_SEGMENT_WBACKSUB   } from '../subworkflows/local/mesmer_segment_wbacksub'
include { MESMER_SEGMENT            } from '../subworkflows/local/mesmer_segment'
include { CELLSAM_SEGMENT_WBACKSUB  } from '../subworkflows/local/cellsam_segment_wbacksub'
include { CELLSAM_SEGMENT           } from '../subworkflows/local/cellsam_segment'
include { SOPA_SEGMENT              } from '../subworkflows/local/sopa_segment'
include { SOPA_SEGMENT_WBACKSUB     } from '../subworkflows/local/sopa_segment_wbacksub'
include { AVITI_SEGMENT             } from '../subworkflows/local/aviti_segment'
include { CELLPOSEMODEL             } from '../modules/local/cellposemodel/main.nf'
include { KRONOS2EMBEDDINGS         } from '../modules/local/kronos2embeddings/main.nf'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_sp_segment_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SP_SEGMENT {

    take:
    ch_samplesheet       // channel: samplesheet read in from --input
    ch_aviti_samplesheet // channel: samplesheet read in from --aviti_input
    main:

    ch_versions = channel.empty()

    //
    // Construct channel for background subtraction/segmentation workflow for MESMER
    //
    ch_samplesheet.branch { it ->
        backsub_only: it[1] == true &&  // run_backsub true
                        it[2] == false && // run_mesmer false
                        it[3] == false && // run_cellpose false
                        it[4] == false    // run_cellsam false
        backsub_mesmer: it[1] == true && it[2] == true // run_backsub true, run_mesmer true
        mesmer_only: it[1] == false && it[2] == true // run_backsub false, run_mesmer true
    }.set { ch_mesmer }

    //
    // Run the BACKGROUNDSUBTRACT subworkflow for samples that ONLY require
    // background subtraction (no segmentation)
    //
    BACKGROUNDSUBTRACT(
        ch_mesmer.backsub_only.map {
            sample,
            _run_backsub,
            _run_mesmer,
            _run_cellpose,
            _run_cellsam,
            tiff,
            _nuclear_channel,
            _membrane_channels -> [
                sample,
                tiff
            ]
        }
    )

    //
    // Run the MESMER_SEGMENT_WBACKSUB subworkflow for samples that require
    // background subtraction and mesmer segmentation
    //
    MESMER_SEGMENT_WBACKSUB(
        ch_mesmer.backsub_mesmer
    )

    //
    // Run MESMER_SEGMENT subworkflow for samples that ONLY require mesmer segmentation
    //
    MESMER_SEGMENT(
        ch_mesmer.mesmer_only
    )

    //
    // Construct channel for only CELLPOSE subworkflow
    //
    ch_samplesheet.filter {
        sample -> sample[3] == true // run_cellpose true for sample
    }.map {
        sample,
        run_backsub,
        _run_mesmer,
        _run_cellpose,
        _run_cellsam,
        tiff,
        nuclear_channel,
        membrane_channels -> [
            sample,
            run_backsub,
            tiff,
            nuclear_channel,
            membrane_channels
        ]
    }.branch { it ->
        with_backsub: it[1] == true // run_backsub true
        no_backsub: it[1] == false  // run_backsub false
    }.set { ch_cellpose_samplesheet }

    //
    // Stage the Cellpose model weights once for the whole run.
    //
    // Invoked ONCE here rather than inside the Cellpose subworkflows, for the
    // same reason KRONOS is hoisted below: SOPA_SEGMENT is instantiated twice
    // (with and without background subtraction) and instantiates
    // SOPA_SEGMENT_COMPARTMENT twice again (nuclear and whole-cell), so a call
    // sited any lower would fetch the weights up to four times and could not
    // guarantee that every patch task in a run read the same file.
    //
    // The trigger channel is the Cellpose samplesheet reduced to the model
    // name: `.first()` emits nothing when no sample runs Cellpose, so Mesmer-
    // and CellSAM-only runs never pay the multi-GB download, and mapping to a
    // constant keeps the `-resume` cache key independent of which sample
    // happened to arrive first. Using the model name as that constant means
    // changing the model correctly invalidates the cached download.
    //
    // A custom model path is loaded directly by cellpose and needs no download,
    // so in that case an empty staged directory is enough.
    //
    if (params.cellpose_models_dir) {
        ch_cellpose_models = channel.value(file(params.cellpose_models_dir, checkIfExists: true))
    } else if (file(params.cellpose_pretrained_model).exists()) {
        ch_cellpose_models = channel.value(file(params.cellpose_pretrained_model).parent)
    } else {
        CELLPOSEMODEL(
            ch_cellpose_samplesheet.with_backsub
                .mix(ch_cellpose_samplesheet.no_backsub)
                .map { params.cellpose_pretrained_model }
                .first()
        )
        // This is a value channel, because the trigger above is one and
        // Nextflow propagates that to the outputs. It has to be: a queue
        // channel would be consumed by the first patch task and starve
        // the rest.
        ch_cellpose_models = CELLPOSEMODEL.out.models
        ch_versions = ch_versions.mix(CELLPOSEMODEL.out.versions)
    }

    //
    // Run CELLPOSE subworkflow for samples that require background subtraction
    //
    SOPA_SEGMENT_WBACKSUB(
        ch_cellpose_samplesheet.with_backsub.map { sample,
            _run_backsub,
            tiff,
            nuclear_channel,
            membrane_channels ->
            [ sample, tiff, nuclear_channel, membrane_channels ]
        },
        ch_cellpose_models
    )

    //
    // Run CELLPOSE subworkflow for samples that ONLY require cellpose segmentation
    //
    SOPA_SEGMENT(
        ch_cellpose_samplesheet.no_backsub.map { sample,
            _run_backsub,
            tiff,
            nuclear_channel,
            membrane_channels ->
            [ sample, tiff, nuclear_channel, membrane_channels ]
        },
        ch_cellpose_models
    )

    //
    // Construct channel for CellSAM segmentation workflow
    //
    ch_samplesheet.branch { it ->
        backsub_cellsam: it[1] == true && it[4] == true // run_backsub true, run_cellsam true
        cellsam_only: it[1] == false && it[4] == true   // run_backsub false, run_cellsam true
    }.set { ch_cellsam }

    //
    // Run the CELLSAM_SEGMENT_WBACKSUB subworkflow for samples that require
    // background subtraction and CellSAM segmentation
    //
    CELLSAM_SEGMENT_WBACKSUB(
        ch_cellsam.backsub_cellsam
    )

    //
    // Run CELLSAM_SEGMENT subworkflow for samples that ONLY require CellSAM segmentation
    //
    CELLSAM_SEGMENT(
        ch_cellsam.cellsam_only
    )

    //
    // Run the AVITI24/Teton/Teton Atlas segmentation subworkflow, if any AVITI
    // runs were supplied via --aviti_input. This is a fully separate path
    // from everything above: it discovers wells/tiles from a run directory
    // rather than taking a pre-segmented sample TIFF, and only its final
    // stitched per-well outputs rejoin the existing CELLMEASUREMENT/KRONOS
    // contract.
    //
    // The custom nuclear model is staged once here (a `path()` input, not a
    // download) for the same reason ch_cellpose_models is hoisted above: it
    // is read by every tile task in the run, and staging inside the
    // subworkflow itself would still be shared correctly, but keeping every
    // "stage a shared model once" decision at this level makes the
    // parallelism/model-sharing shape of the whole workflow easier to read
    // from one place.
    //
    ch_aviti_nuclear_model = params.aviti_nuclear_model_path
        ? channel.value(file(params.aviti_nuclear_model_path, checkIfExists: true))
        : channel.empty()

    AVITI_SEGMENT(
        ch_aviti_samplesheet,
        ch_aviti_nuclear_model
    )
    ch_versions = ch_versions.mix(AVITI_SEGMENT.out.versions)

    //
    // Optional KRONOS embedding extraction
    //
    // Invoked ONCE here rather than inside each segmentation subworkflow.
    // KRONOS is segmenter-agnostic, and a cohort-wide pass must see every
    // sample: a channel operator inside e.g. MESMER_SEGMENT would only ever
    // observe the Mesmer samples, so a cohort-level statistic computed there
    // would silently differ per segmenter.
    //
    if (params.enable_kronos) {

        ch_kronos_input = MESMER_SEGMENT.out.kronos_input
            .mix(MESMER_SEGMENT_WBACKSUB.out.kronos_input)
            .mix(SOPA_SEGMENT.out.kronos_input)
            .mix(SOPA_SEGMENT_WBACKSUB.out.kronos_input)
            .mix(CELLSAM_SEGMENT.out.kronos_input)
            .mix(CELLSAM_SEGMENT_WBACKSUB.out.kronos_input)
            .mix(AVITI_SEGMENT.out.kronos_input)

        ch_kronos_annotations = MESMER_SEGMENT.out.annotations
            .mix(MESMER_SEGMENT_WBACKSUB.out.annotations)
            .mix(SOPA_SEGMENT.out.annotations)
            .mix(SOPA_SEGMENT_WBACKSUB.out.annotations)
            .mix(CELLSAM_SEGMENT.out.annotations)
            .mix(CELLSAM_SEGMENT_WBACKSUB.out.annotations)
            .mix(AVITI_SEGMENT.out.annotations)

        // The samplesheet's per-sample nuclear channel is the model's DAPI hint
        // (preferred_dapi), which is KRONOS2's only marker-alias mechanism.
        ch_nuclear = ch_samplesheet.map {
            meta,
            _run_backsub,
            _run_mesmer,
            _run_cellpose,
            _run_cellsam,
            _tiff,
            nuclear_channel,
            _membrane_channels -> [ meta, nuclear_channel ]
        }
        // AVITI's nuclear channel is always named "Nucleus" -- fixed by the
        // instrument's own cell-paint channel naming, unlike the free-form
        // per-sample mapping the COMET/MIBI samplesheet supplies.
        .mix(
            AVITI_SEGMENT.out.kronos_input.map { meta, _tiff, _whole_cell_mask -> [ meta, 'Nucleus' ] }
        )

        // Everything is joined on meta before the call: Nextflow consumes
        // separate channels positionally, so mixing six upstream sources that
        // finish at different times could otherwise pair one sample's image
        // with another sample's annotations.
        //
        // KRONOS2 derives cells from the GeoJSON polygons, so the whole-cell
        // mask carried by kronos_input is not needed here and is dropped
        // before staging.
        KRONOS2EMBEDDINGS(
            ch_kronos_input
                .join(ch_kronos_annotations, by: 0)
                .join(ch_nuclear, by: 0)
                .map { meta, tiff, _whole_cell_mask, geojson, nuclear_channel ->
                    [ meta, tiff, geojson, nuclear_channel ]
                },
            file(params.kronos_model_path)
        )
        ch_versions = ch_versions.mix(KRONOS2EMBEDDINGS.out.versions.first())
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  ''  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_collated_versions     // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
