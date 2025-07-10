/*
 * This subworkflow uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */

process SOPACONVERT {
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.0.7'
        : 'docker.io/quentinblampey/sopa:2.0.7'}"

    input:
    tuple val(meta), path(tiff)

    output:
    tuple val(meta), path("*.zarr"), emit: spatial_data

    script:
    def args = task.ext.args ?: ''
    """
    sopa convert \\
        ${args} \\
        --sdata-path ${meta.id}.zarr \\
        --technology ${params.technology} \\
        ${tiff}
    """
}

workflow SOPASEGMENT {

    take:
    ch_sopa // channel: [ (meta, tiff, nuclear_channel, membrane_channels) ]

    main:

    ch_versions = Channel.empty()

    ch_sopa.map { meta, tiff, nuclear_channel, membrane_channels, skip_measurements ->
        [ meta, tiff ]
    }.set { ch_convert }

    //
    // Run SOPA convert to convert tiff to zarr format
    //
    SOPACONVERT(
        ch_convert
    )

    emit:
    zarr     = SOPACONVERT.out.spatial_data    // channel: [ val(meta), *.zarr ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
