process MESMERSEGMENT {
    tag "$meta.id"
    label 'process_multi'
    secret 'DEEPCELL_ACCESS_TOKEN'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/wehi-soda-hub/mesmersegmentation:0.2.0'

    input:
    tuple val(meta), path(tiff), val(nuclear_channel), val(membrane_channels)
    val(compartment)

    output:
    tuple val(meta), path("*.tiff"), emit: segmentation_mask
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def membrane_channel_args = membrane_channels != '' && membrane_channels != [] ?
        membrane_channels.split(":").collect {
            "--membrane-channel \"${it}\""
        }.join(' ') : ''
    """
    mesmer-segment \\
        ${tiff} \\
        --compartment ${compartment} \\
        --nuclear-channel ${nuclear_channel} \\
        ${membrane_channel_args} \\
        ${args} \\
        > "${prefix}_${compartment}.tiff"

    # Post-process: transpose mask if dimensions are swapped
    mesmer_postprocess.py "${tiff}" "${prefix}_${compartment}.tiff" ${params.mesmer_transpose_mask ? '--force-transpose' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mesmersegmentation: v0.1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.tiff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mesmersegmentation: v0.1.0
    END_VERSIONS
    """
}
