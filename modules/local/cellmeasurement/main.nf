process CELLMEASUREMENT {
    tag "$meta.id"
    label 'process_medium_memory'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/wehi-soda-hub/cellmeasurement:latest"

    input:
    tuple val(meta),
        path(tiff),
        path(nuclear_mask),
        path(whole_cell_mask),
        val(cell_expansion),
        val(skip_measurements)
    val(pixel_size_microns)

    output:
    tuple val(meta), path("*.geojson"), emit: annotations
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def cell_expansion_arg = cell_expansion.first() != [] ? "--cell-expansion ${cell_expansion.first()}" : ''
    def skip_measurements_arg = skip_measurements ? '--skip-measurements=true' : '--skip-measurements=false'
    """
    /entrypoint.sh \\
        --args="${args} \\
            --nuclear-mask=${nuclear_mask} \\
            --whole-cell-mask=${whole_cell_mask} \\
            --tiff-file=${tiff} \\
            --output-file=${prefix}.geojson \\
            ${cell_expansion_arg} \\
            ${skip_measurements_arg}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellmeasurement: \$(/entrypoint.sh --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.geojson

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellmeasurement: \$(/entrypoint.sh --version)
    END_VERSIONS
    """
}
