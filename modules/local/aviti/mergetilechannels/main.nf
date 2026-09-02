process AVITIMERGETILECHANNELS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Same tifffile/numpy/typer image already used by PARQUETTOTIFF.
    container 'community.wave.seqera.io/library/tifffile_pyarrow_rasterio_shapely_pruned:9cab11ac36e81144'

    input:
    tuple val(meta), path(nucleus_tif), path(membrane_tif), path(actin_tif)

    output:
    tuple val(meta), path("${meta.id}_image.ome.tif"), emit: image
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // actin_tif is staged as a literal "NO_FILE" placeholder in 2-channel mode
    // (see the subworkflow), since Nextflow tuples cannot carry an optional
    // path slot directly.
    def actin_arg = (actin_tif.name != 'NO_FILE') ? "--actin-tif ${actin_tif}" : ''
    """
    aviti_merge_tile_channels.py \\
        --nucleus-tif ${nucleus_tif} \\
        --membrane-tif ${membrane_tif} \\
        ${actin_arg} \\
        --output ${meta.id}_image.ome.tif \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        tifffile: \$(python3 -c "import tifffile; print(tifffile.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_image.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        tifffile: \$(python3 -c "import tifffile; print(tifffile.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
