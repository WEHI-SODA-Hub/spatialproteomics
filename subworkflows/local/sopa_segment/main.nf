include { COMBINECHANNELS                                    } from '../../../modules/local/combinechannels/main.nf'
include { SOPA_SEGMENT_COMPARTMENT as SOPA_SEGMENT_NUCLEAR   } from '../sopa_segment_compartment/main.nf'
include { SOPA_SEGMENT_COMPARTMENT as SOPA_SEGMENT_WHOLECELL } from '../sopa_segment_compartment/main.nf'
include { CELLMEASUREMENT                                    } from '../../../modules/local/cellmeasurement/main.nf'
include { KRONOSEMBEDDINGS                                    } from '../../../modules/local/kronosembeddings/main.nf'
include { SEGMENTATIONREPORT                                 } from '../../../modules/local/segmentationreport/main.nf'

workflow SOPA_SEGMENT {

    take:
    ch_sopa // channel: [ (meta, tiff, nuclear_channel, membrane_channels) ]

    main:

    ch_versions = Channel.empty()

    //
    // Combine membrane channels into a single channel
    //
    COMBINECHANNELS(
        ch_sopa
    )
    ch_versions = ch_versions.mix(COMBINECHANNELS.out.versions.first())

    // Replace tiff with combined_tiff
    // If there are multiple membrane channels, rename to 'combined_membrane'
    COMBINECHANNELS.out.combined_tiff
        .join( ch_sopa, by: 0 )
        .map { meta, combined_tiff, _tiff, nuclear_channel, membrane_channels ->
            def membrane_name = membrane_channels.split(':').size() == 1 ?
                membrane_channels : 'combined_membrane'
            [ meta, combined_tiff, nuclear_channel, membrane_name ]
        }.set { ch_combined }

    //
    // Run segmentation for nuclear compartment
    //
    SOPA_SEGMENT_NUCLEAR(
        ch_combined.map {
            meta,
            tiff,
            nuclear_channel,
            _membrane_channels -> [
                meta,
                tiff,
                nuclear_channel,
                '' // no membrane channels for nuclear segmentation
            ]
        },
        'nuclear'
    )
    ch_versions = ch_versions.mix(SOPA_SEGMENT_NUCLEAR.out.versions.first())

    //
    // Run segmentation for whole-cell compartment
    //
    SOPA_SEGMENT_WHOLECELL(
        ch_combined,
        'whole-cell'
    )
    ch_versions = ch_versions.mix(SOPA_SEGMENT_WHOLECELL.out.versions.first())

    //
    // Create a channel for cell measurement
    //
    SOPA_SEGMENT_NUCLEAR.out.tiff
        .join(SOPA_SEGMENT_WHOLECELL.out.tiff, by: 0)
        .join(ch_sopa, by: 0)
        .map {
            meta,
            nuclear_tiff,
            wholecell_tiff,
            tiff,
            _nuc_chan,
            _mem_chans -> [
                meta,
                tiff,
                nuclear_tiff,
                wholecell_tiff
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

        // Create channel for KRONOS input: original tiff + whole-cell mask
        ch_sopa
            .join(SOPA_SEGMENT_WHOLECELL.out.tiff)
            .map {
                meta,
                tiff,
                _nuclear_channel,
                _membrane_channels,
                whole_cell_mask -> [
                    meta,
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
        ch_combined
            .join(CELLMEASUREMENT.out.annotations)
            .map {
                sample,
                combined_tiff,
                nuclear_channel,
                membrane_channels,
                annotations -> [
                    sample,
                    annotations,
                    false, // run_mesmer
                    true,  // run_cellpose
                    nuclear_channel.first(),
                    membrane_channels.first(),
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

        ch_report = SEGMENTATIONREPORT.out.report  // channel: [ val(meta), *.html ]
    }

    emit:
    annotations          = CELLMEASUREMENT.out.annotations   // channel: [ val(meta), *.geojson ]
    kronos_embeddings    = ch_kronos_embeddings               // channel: [ val(meta), *.csv ] OPTIONAL
    kronos_marker_report = ch_kronos_marker_report            // channel: [ val(meta), *.txt ] OPTIONAL
    kronos_merged_geojson = ch_kronos_merged_geojson           // channel: [ val(meta), *.geojson ] OPTIONAL
    report               = ch_report                         // channel: [ val(meta), *.html ] OPTIONAL

    versions = ch_versions                          // channel: [ versions.yml ]
}
