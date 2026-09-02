process AVITINUCLEARSEGMENT {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://community.wave.seqera.io/library/python_numpy_pandas_scikit-image_pruned:4cec16761766b006'
        : 'community.wave.seqera.io/library/python_numpy_pandas_scikit-image_pruned:4cec16761766b006'}"

    input:
    tuple val(meta), path(nucleus_tif)
    path model_path

    output:
    // Named after the source tile, not meta.id: this is exactly the AVITI
    // viewer filename (`<tile>_Nuclear.tif`) expected under `Well<well>/`.
    tuple val(meta), path("${meta.tile}_Nuclear.tif"), emit: nuclear_mask
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    aviti_nuclear_segment.py \\
        ${nucleus_tif} \\
        --output ${meta.tile}_Nuclear.tif \\
        --model-path ${model_path} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(python3 -c "from importlib.metadata import version; print(version('cellpose'))")
        torch: \$(python3 -c "import torch; print(torch.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tile}_Nuclear.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(python3 -c "from importlib.metadata import version; print(version('cellpose'))")
        torch: \$(python3 -c "import torch; print(torch.__version__)")
    END_VERSIONS
    """
}
