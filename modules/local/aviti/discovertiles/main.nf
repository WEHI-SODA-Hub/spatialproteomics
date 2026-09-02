process AVITIDISCOVERTILES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Reuses the tifffile/pyarrow/typer image already pulled for
    // PARQUETTOTIFF: this module's own runtime needs is just typer (json/csv
    // are stdlib), and that container already satisfies it without pulling a
    // new one.
    container 'community.wave.seqera.io/library/tifffile_pyarrow_rasterio_shapely_pruned:9cab11ac36e81144'

    input:
    tuple val(meta), path(run_dir)

    output:
    tuple val(meta), path("${meta.id}.manifest.csv"), emit: manifest
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // meta.wells follows this pipeline's colon-separated list convention
    // (matching membrane_channels in the COMET/MIBI samplesheet), but the
    // underlying script takes a comma-separated list -- translate here
    // rather than adding a second list separator to the script's own CLI.
    def wells_arg = meta.wells ? "--wells '${meta.wells.replace(':', ',')}'" : ''
    """
    aviti_discover_tiles.py \\
        ${run_dir} \\
        --output ${meta.id}.manifest.csv \\
        ${wells_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        typer: \$(python3 -c "import typer; print(typer.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo "well,tile,x_mm,y_mm,nucleus_tif,membrane_tif,actin_tif,channel_mode" > ${meta.id}.manifest.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        typer: \$(python3 -c "import typer; print(typer.__version__)")
    END_VERSIONS
    """
}
