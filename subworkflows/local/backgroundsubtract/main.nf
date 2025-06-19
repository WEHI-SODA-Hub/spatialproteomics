include { EXTRACTMARKERS      } from '../../../modules/local/extractmarkers/main.nf'
include { BACKSUB             } from '../../../modules/nf-core/backsub/main.nf'

workflow BACKGROUNDSUBTRACT {

    take:
    ch_backsub // channel: background subtraction parameters (sample name and tiff)

    main:

    ch_versions = Channel.empty()

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

    emit:
    backsub_tif   = BACKSUB.out.backsub_tif    // channel: [ val(meta), *.ome.tif ]
    markers       = BACKSUB.out.markerout      // channel: [ val(meta2), markers.csv ]

    versions      = ch_versions                     // channel: [ versions.yml ]
}
