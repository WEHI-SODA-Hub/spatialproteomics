process KRONOS2EMBEDDINGS {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    // TODO: build via Seqera Wave from environment.yml and pin the digest here,
    // as was done for the KRONOS1 module. The cu124 torch/xformers stack has no
    // suitable public image.
    container 'community.wave.seqera.io/library/kronos2embeddings:PLACEHOLDER'

    input:
    tuple val(meta), path(tiff), path(geojson, stageAs: 'cellmeas_input/*'), val(nuclear_channel)
    path(kronos2_model)

    output:
    tuple val(meta), path("*.geojson{,.gz}")     , emit: annotations
    tuple val(meta), path("*_marker_report.txt") , emit: marker_report
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def geojson_ext = geojson.name.endsWith('.gz') ? '.geojson.gz' : '.geojson'
    // preferred_dapi is the model's only alias mechanism, so a slide stained
    // with DRAQ5/Hoechst needs it to receive DAPI stats rather than defaults.
    // The global param overrides the samplesheet's per-sample nuclear channel.
    def nuclear = params.kronos_nuclear_marker ?: nuclear_channel
    def nuclear_arg = nuclear ? "--nuclear-marker '${nuclear}'" : ''
    """
    kronos2_embeddings.py \\
        --tiff ${tiff} \\
        --geojson ${geojson} \\
        --model-path ${kronos2_model} \\
        --output-geojson ${prefix}${geojson_ext} \\
        --marker-report ${prefix}_marker_report.txt \\
        ${nuclear_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pytorch: \$(python -c "import torch; print(torch.__version__)")
        transformers: \$(python -c "import transformers; print(transformers.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def geojson_ext = geojson.name.endsWith('.gz') ? '.geojson.gz' : '.geojson'
    """
    touch ${prefix}${geojson_ext}
    touch ${prefix}_marker_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pytorch: \$(python -c "import torch; print(torch.__version__)")
        transformers: \$(python -c "import transformers; print(transformers.__version__)")
    END_VERSIONS
    """
}
