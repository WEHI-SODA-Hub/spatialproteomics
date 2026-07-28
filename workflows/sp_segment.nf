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
include { KRONOSEMBEDDINGS          } from '../modules/local/kronosembeddings/main.nf'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_sp_segment_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SP_SEGMENT {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

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
        it[3] == true // run_cellpose true for sample
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
        with_backsub: it[1] == true// run_backsub true
        no_backsub: it[1] == false // run_backsub false
    }.set { ch_cellpose_samplesheet }

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
        }
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
        }
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

        ch_kronos_annotations = MESMER_SEGMENT.out.annotations
            .mix(MESMER_SEGMENT_WBACKSUB.out.annotations)
            .mix(SOPA_SEGMENT.out.annotations)
            .mix(SOPA_SEGMENT_WBACKSUB.out.annotations)
            .mix(CELLSAM_SEGMENT.out.annotations)
            .mix(CELLSAM_SEGMENT_WBACKSUB.out.annotations)

        // KRONOSEMBEDDINGS takes the image/mask and the GeoJSON as two separate
        // channels, which Nextflow consumes positionally. Join on meta first and
        // derive both inputs from that single joined channel, so a sample's mask
        // can never be paired with another sample's annotations.
        ch_kronos_joined = ch_kronos_input.join(ch_kronos_annotations, by: 0)

        KRONOSEMBEDDINGS(
            ch_kronos_joined.map { meta, tiff, whole_cell_mask, _geojson ->
                [ meta, tiff, whole_cell_mask ]
            },
            file(params.kronos_model_path),
            file(params.kronos_marker_metadata),
            ch_kronos_joined.map { meta, _tiff, _whole_cell_mask, geojson ->
                [ meta, geojson ]
            }
        )
        ch_versions = ch_versions.mix(KRONOSEMBEDDINGS.out.versions.first())
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
