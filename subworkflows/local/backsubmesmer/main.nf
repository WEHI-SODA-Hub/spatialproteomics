include { EXTRACTMARKERS      } from '../../../modules/local/extractmarkers/main.nf'
include { BACKSUB             } from '../../../modules/nf-core/backsub/main.nf'

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
    // Run background subtraction module on tiff with extracted markers
    //
    BACKSUB(
        ch_backsub,
        EXTRACTMARKERS.out.markers
    )

    // TODO: run mesmer module here

    emit:
    backsub_tif   = BACKSUB.out.backsub_tif    // channel: [ val(meta), *.ome.tif ]
    markers       = BACKSUB.out.markerout      // channel: [ val(meta2), markers.csv ]

    versions      = ch_versions                     // channel: [ versions.yml ]
}
