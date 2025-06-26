include { EXTRACTMARKERS             } from '../../../modules/local/extractmarkers/main.nf'
include { BACKSUB                    } from '../../../modules/nf-core/backsub/main.nf'
include { MESMERSEGMENT as MESMERWC  } from '../../../modules/local/mesmersegment/main.nf'
include { MESMERSEGMENT as MESMERNUC } from '../../../modules/local/mesmersegment/main.nf'

workflow BACKSUBMESMER {

    take:
    ch_backsub_mesmer // channel: segmentation parameters

    main:

    ch_versions = Channel.empty()

    ch_backsub_mesmer.map {
        sample,
        run_backsub,
        run_mesmer,
        seg_tiff,
        seg_nuclear_channel,
        seg_membrane_channels,
        seg_combine_method,
        seg_level,
        seg_maxima_threshold,
        seg_interior_threshold,
        seg_maxima_smooth,
        seg_min_nuclei_area,
        seg_remove_border_cells,
        seg_pixel_expansion,
        seg_padding -> [
            sample,
            seg_tiff,
        ]
    }.set { ch_backsub }

    //
    // Extract markers from the input tiff file
    //
    EXTRACTMARKERS(
        ch_backsub
    )

    //
    // Run BACKSUB module on tiff with extracted markers
    //
    BACKSUB(
        ch_backsub,
        EXTRACTMARKERS.out.markers
    )

    // Prepare channel for input to mesmer
    BACKSUB.out.backsub_tif
        .join(ch_backsub_mesmer)

    ch_backsub_mesmer
        .join(BACKSUB.out.backsub_tif)
        .map {
            sample,
            run_backsub,
            run_mesmer,
            seg_tiff,
            seg_nuclear_channel,
            seg_membrane_channels,
            seg_combine_method,
            seg_level,
            seg_maxima_threshold,
            seg_interior_threshold,
            seg_maxima_smooth,
            seg_min_nuclei_area,
            seg_remove_border_cells,
            seg_pixel_expansion,
            seg_padding,
            backsub_tif -> [
                sample,
                backsub_tif,
                seg_nuclear_channel,
                seg_membrane_channels,
                seg_combine_method,
                seg_level,
                seg_maxima_threshold,
                seg_interior_threshold,
                seg_maxima_smooth,
                seg_min_nuclei_area,
                seg_remove_border_cells,
                seg_pixel_expansion,
                seg_padding
            ]
        }.set { ch_mesmer }


    ch_mesmer.view()
    //
    // Run MESMERSEGMENT module on the background subtracted tiff
    // for whole-cell segmentation
    //
    MESMERWC(
        ch_mesmer,
        "whole-cell"
    )

    //
    // Run MESMERSEGMENT module as above, but this time for nuclear segmentation
    //
    MESMERNUC(
        ch_mesmer,
        "nuclear"
    )


    emit:
    whole_cell_tif   = MESMERWC.out.segmentation_mask    // channel: [ val(meta), *.tiff ]
    nuclear_tif      = MESMERNUC.out.segmentation_mask   // channel: [ val(meta), *.tiff ]

    versions         = ch_versions                     // channel: [ versions.yml ]
}
