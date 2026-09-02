process AVITISTITCHWELL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // Same tifffile/numpy/typer image already used by PARQUETTOTIFF.
    container 'community.wave.seqera.io/library/tifffile_pyarrow_rasterio_shapely_pruned:9cab11ac36e81144'

    input:
    // tile_rows is a list of maps ([tile, x_mm, y_mm, cell_mask, nuclear_mask,
    // image_tif]) built by the subworkflow; cell_masks/nuclear_masks/images
    // are the corresponding staged files, named exactly as referenced in
    // tile_rows so the manifest built below resolves them.
    tuple val(meta), val(tile_rows), path(cell_masks), path(nuclear_masks), path(images)

    output:
    tuple val(meta), path("${meta.id}_cell_stitched.tif")   , emit: cell_mask
    tuple val(meta), path("${meta.id}_nuclear_stitched.tif"), emit: nuclear_mask
    tuple val(meta), path("${meta.id}_image_stitched.tif")  , emit: image
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def manifest_lines = (
        ['tile,x_mm,y_mm,cell_mask,nuclear_mask,image_tif'] +
        tile_rows.collect { r -> "${r.tile},${r.x_mm},${r.y_mm},${r.cell_mask},${r.nuclear_mask},${r.image_tif}" }
    ).join('\n    ')
    """
    cat > manifest.csv <<'AVITI_MANIFEST_EOF'
    ${manifest_lines}
    AVITI_MANIFEST_EOF

    aviti_stitch.py \\
        manifest.csv \\
        --output-cell ${meta.id}_cell_stitched.tif \\
        --output-nuclear ${meta.id}_nuclear_stitched.tif \\
        --output-image ${meta.id}_image_stitched.tif \\
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
    touch ${meta.id}_cell_stitched.tif
    touch ${meta.id}_nuclear_stitched.tif
    touch ${meta.id}_image_stitched.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        tifffile: \$(python3 -c "import tifffile; print(tifffile.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
