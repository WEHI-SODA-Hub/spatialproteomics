process KRONOSEMBEDDINGS {
    tag "$meta.id"
    label 'process_multi'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/python_git_numpy_pandas_pruned:16c3de943194d1fd'

    input:
    tuple val(meta),
        path(tiff),
        path(whole_cell_mask)
    path(kronos_model)
    path(marker_metadata)
    tuple val(meta2), path(geojson)

    output:
    tuple val(meta), path("*_kronos_embeddings.csv")  , emit: embeddings
    tuple val(meta), path("*_marker_report.txt")       , emit: marker_report
    tuple val(meta), path("*_kronos_merged.geojson")   , emit: merged_geojson, optional: true
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def merge_args = params.kronos_merge_geojson ? "--geojson ${geojson} --merge-geojson" : ''
    """
    kronos_embeddings.py \\
        --tiff ${tiff} \\
        --mask ${whole_cell_mask} \\
        --model-path ${kronos_model}/*.pt \\
        --marker-metadata ${marker_metadata} \\
        --output ${prefix}_kronos_embeddings.csv \\
        --sample-id ${prefix} \\
        ${merge_args} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pytorch: \$(python -c "import torch; print(torch.__version__)")
        kronos: \$(python -c "import kronos; print('0.1.0')")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_kronos_embeddings.csv
    touch ${prefix}_kronos_embeddings_marker_report.txt
    ${params.kronos_merge_geojson ? "touch ${prefix}_kronos_merged.geojson" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pytorch: \$(python -c "import torch; print(torch.__version__)")
        kronos: \$(python -c "import kronos; print('0.1.0')")
    END_VERSIONS
    """
}
