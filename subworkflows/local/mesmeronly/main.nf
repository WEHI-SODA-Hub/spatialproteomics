include { MESMERSEGMENT as MESMERWC  } from '../../../modules/local/mesmersegment/main.nf'
include { MESMERSEGMENT as MESMERNUC } from '../../../modules/local/mesmersegment/main.nf'
include { CELLMEASUREMENT            } from '../../../modules/local/cellmeasurement/main.nf'

workflow MESMERONLY {

    take:
    ch_mesmeronly // channel: segmentation parameters

    main:

    ch_versions = Channel.empty()

    ch_mesmeronly.map {
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
        seg_skip_measurements -> [
            sample,
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
            seg_padding
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

    //
    // Run MESMERSEGMENT module as above, but this time for nuclear segmentation
    //
    MESMERNUC(
        ch_mesmer,
        "nuclear"
    )

    // Create channel for CELLMEASUREMENT input adding the segmentation masks
    ch_mesmeronly
        .join(MESMERNUC.out.segmentation_mask)
        .join(MESMERWC.out.segmentation_mask)
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
            seg_skip_measurements,
            nuclear_mask,
            whole_cell_mask -> [
                sample,
                seg_tiff,
                nuclear_mask,
                whole_cell_mask,
                seg_skip_measurements
            ]
        }.set { ch_cellmeasurement }

    //
    // Run CELLMEASUREMENT module on the whole-cell and nuclear segmentation masks
    //
    CELLMEASUREMENT(
        ch_cellmeasurement,
        params.pixel_size_microns
    )

    emit:
    annotations      = CELLMEASUREMENT.out.annotations   // channel: [ val(meta), *.geojson ]
    whole_cell_tif   = MESMERWC.out.segmentation_mask    // channel: [ val(meta), *.tiff ]
    nuclear_tif      = MESMERNUC.out.segmentation_mask   // channel: [ val(meta), *.tiff ]

    versions         = ch_versions                       // channel: [ versions.yml ]
}
