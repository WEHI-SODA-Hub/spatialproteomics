process MESMERSEGMENT {
    tag "$meta.id"
    label 'process_medium_memory'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/wehi-soda-hub/mibisegmentation:latest'

    input:
    tuple val(meta),
        path(tiff),
        val(nuclear_channel),
        val(membrane_channels),
        val(combine_method),
        val(level),
        val(maxima_threshold),
        val(interior_threshold),
        val(maxima_smooth),
        val(min_nuclei_area),
        val(remove_border_cells),
        val(pixel_expansion),
        val(padding)
    val(compartment)

    output:
    tuple val(meta), path("*.tiff"), emit: segmentation_mask
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def membrane_channel_args = membrane_channels.first().split(":")
        .collect { "--membrane-channel ${it}" }.join(' ')

    // TODO: These params are a mess now, consider refactoring
    def combine_method_arg = combine_method.first() != [] ? "--combine-method ${combine_method.first()}" : ''
    def level_arg = level.first() != [] ? "--segmentation-level ${level.first()}" : ''
    def maxima_threshold_arg = maxima_threshold.first() != [] ? "--maxima-threshold ${maxima_threshold.first()}" : ''
    def interior_threshold_arg = interior_threshold.first() != [] ? "--interior-threshold ${interior_threshold.first()}" : ''
    def maxima_smooth_arg = maxima_smooth.first() != [] ? "--maxima-smooth ${maxima_smooth.first()}" : ''
    def min_nuclei_area_arg = min_nuclei_area.first() != [] ? "--min-nuclei-area ${min_nuclei_area.first()}" : ''
    def remove_border_cells_arg = remove_border_cells.first() ? '--remove-cells-touching-border' : '--no-remove-cells-touching-border'
    def pixel_expansion_arg = pixel_expansion.first() != [] ? "--pixel-expansion ${pixel_expansion.first()}" : ''
    def padding_arg = padding.first() != [] ? "--padding ${padding.first()}" : ''
    """
    mesmer-segment \\
        ${tiff} \\
        --compartment ${compartment} \\
        --nuclear-channel ${nuclear_channel.first()} \\
        ${membrane_channel_args} \\
        ${combine_method_arg} \\
        ${level_arg} \\
        ${maxima_threshold_arg} \\
        ${interior_threshold_arg} \\
        ${maxima_smooth_arg} \\
        ${min_nuclei_area_arg} \\
        ${remove_border_cells_arg} \\
        ${pixel_expansion_arg} \\
        ${padding_arg} \\
        $args \\
        > "${prefix}_${compartment}.tiff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mesmersegmentation: v0.1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.tiff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mesmersegmentation: v0.1.0
    END_VERSIONS
    """
}
