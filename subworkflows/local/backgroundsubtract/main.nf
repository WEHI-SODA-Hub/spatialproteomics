include { EXTRACTMARKERS      } from '../../../modules/local/extractmarkers/main.nf'

workflow BACKGROUNDSUBTRACT {

    take:
    ch_backsub // channel: background subtraction parameters (sample name and tiff)

    main:

    ch_versions = Channel.empty()

    EXTRACTMARKERS(
        ch_backsub
    )

    emit:
    markers      = EXTRACTMARKERS.out.markers  // channel: [ val(meta), markers.csv ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
