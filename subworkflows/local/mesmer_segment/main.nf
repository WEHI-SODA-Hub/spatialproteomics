include { MESMERSEGMENT as MESMERWC  } from '../../../modules/local/mesmersegment/main.nf'
include { MESMERSEGMENT as MESMERNUC } from '../../../modules/local/mesmersegment/main.nf'
include { CELLMEASUREMENT            } from '../../../modules/local/cellmeasurement/main.nf'
include { KRONOSEMBEDDINGS            } from '../../../modules/local/kronosembeddings/main.nf'
include { COMBINECHANNELS            } from '../../../modules/local/combinechannels/main.nf'
include { SEGMENTATIONREPORT         } from '../../../modules/local/segmentationreport/main.nf'

workflow MESMER_SEGMENT {

    take:
    ch_mesmer_segment // channel: [ (sample, run_backsub, run_mesmer, run_cellpose, run_cellsam, tiff, nuclear_channel, membrane_channels) ]

    main:

    ch_versions = Channel.empty()

    ch_mesmer_segment.map {
        sample,
        _run_backsub,
        _run_mesmer,
        _run_cellpose,
        _run_cellsam,
        tiff,
        nuclear_channel,
        membrane_channels -> [
            sample,
            tiff,
            nuclear_channel,
            membrane_channels
        ]
    }.set { ch_mesmer }


    //
    // Run MESMERSEGMENT module on the background subtracted tiff
    // for whole-cell segmentation
    //
    MESMERWC(
        ch_mesmer,
        "whole-cell"
    )
    ch_versions = ch_versions.mix(MESMERWC.out.versions.first())


    //
    // Run MESMERSEGMENT module as above, but this time for nuclear segmentation
    //
    MESMERNUC(
        ch_mesmer,
        "nuclear"
    )
    ch_versions = ch_versions.mix(MESMERNUC.out.versions.first())

    // Create channel for CELLMEASUREMENT input adding the segmentation masks
    ch_mesmer_segment
        .join(MESMERNUC.out.segmentation_mask)
        .join(MESMERWC.out.segmentation_mask)
        .map {
            sample,
            _run_backsub,
            _run_mesmer,
            _run_cellpose,
            _run_cellsam,
            tiff,
            _nuclear_channel,
            _membrane_channels,
            nuclear_mask,
            whole_cell_mask -> [
                sample,
                tiff,
                nuclear_mask,
                whole_cell_mask
            ]
        }.set { ch_cellmeasurement }

    //
    // Run CELLMEASUREMENT module on the whole-cell and nuclear segmentation masks
    //
    CELLMEASUREMENT(
        ch_cellmeasurement
    )
    ch_versions = ch_versions.mix(CELLMEASUREMENT.out.versions.first())

    //
    // Optional KRONOS embedding extraction
    //
    ch_kronos_embeddings = Channel.empty()
    ch_kronos_marker_report = Channel.empty()
    ch_kronos_merged_geojson = Channel.empty()
    if (!params.skip_kronos) {

        // Create channel for KRONOS input: tiff + whole-cell mask + geojson
        ch_mesmer_segment
            .join(MESMERWC.out.segmentation_mask)
            .map {
                sample,
                _run_backsub,
                _run_mesmer,
                _run_cellpose,
                _run_cellsam,
                tiff,
                _nuclear_channel,
                _membrane_channels,
                whole_cell_mask -> [
                    sample,
                    tiff,
                    whole_cell_mask
                ]
            }.set { ch_kronos_input }

        KRONOSEMBEDDINGS(
            ch_kronos_input,
            file(params.kronos_model_path),
            file(params.kronos_marker_metadata),
            CELLMEASUREMENT.out.annotations
        )
        ch_versions = ch_versions.mix(KRONOSEMBEDDINGS.out.versions.first())
        ch_kronos_embeddings = KRONOSEMBEDDINGS.out.embeddings
        ch_kronos_marker_report = KRONOSEMBEDDINGS.out.marker_report
        ch_kronos_merged_geojson = KRONOSEMBEDDINGS.out.merged_geojson
    }

    // Optional SEGMENTATIONREPORT module
    ch_report = Channel.empty()
    if (params.generate_report) {

        //
        // Combine channels for report background image
        //
        COMBINECHANNELS(
            ch_mesmer
        )
        ch_versions = ch_versions.mix(COMBINECHANNELS.out.versions.first())

        ch_mesmer_segment
            .join(CELLMEASUREMENT.out.annotations)
            .join(COMBINECHANNELS.out.combined_tiff, by: 0)
            .map {
                sample,
                _run_backsub,
                run_mesmer,
                run_cellpose,
                _run_cellsam,
                _tiff,
                nuclear_channel,
                membrane_channels,
                annotations,
                combined_tiff -> [
                    sample,
                    annotations,
                    run_mesmer,
                    run_cellpose,
                    nuclear_channel,
                    membrane_channels,
                    combined_tiff
                ]
            }.set { ch_segmentationreport }

        //
        // Run SEGMENTATIONREPORT module to generate a report of the segmentation results
        //
        SEGMENTATIONREPORT(
            ch_segmentationreport
        )
        ch_versions = ch_versions.mix(SEGMENTATIONREPORT.out.versions.first())

        ch_report = SEGMENTATIONREPORT.out.report
    }

    emit:
    annotations      = CELLMEASUREMENT.out.annotations   // channel: [ val(meta), *.geojson ]
    whole_cell_tif   = MESMERWC.out.segmentation_mask    // channel: [ val(meta), *.tiff ]
    nuclear_tif      = MESMERNUC.out.segmentation_mask   // channel: [ val(meta), *.tiff ]
    kronos_embeddings     = ch_kronos_embeddings          // channel: [ val(meta), *.csv ] OPTIONAL
    kronos_marker_report  = ch_kronos_marker_report       // channel: [ val(meta), *.txt ] OPTIONAL
    kronos_merged_geojson = ch_kronos_merged_geojson       // channel: [ val(meta), *.geojson ] OPTIONAL
    report           = ch_report                         // channel: [ val(meta), *.html ] OPTIONAL

    versions         = ch_versions                       // channel: [ versions.yml ]
}
