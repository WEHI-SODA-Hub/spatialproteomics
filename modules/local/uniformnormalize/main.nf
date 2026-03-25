process UNIFORMNORMALIZE {
    tag "uniform"
    label 'process_multi'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/python_tifffile_scikit-image_scikit-learn_pruned:593e00ba324c12b3'

    input:
    path(geojsons)

    output:
    path "*_uniform.geojson"            , emit: normalized_annotations
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def geojson_list = geojsons.join(' ')
    """
    python uniform_normalize_geojson.py \\
        --inputs ${geojson_list} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \\$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    for f in ${geojsons.join(' ')}; do
        base=\\$(basename "\\$f" .geojson)
        touch "\\${base}_uniform.geojson"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \\$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
