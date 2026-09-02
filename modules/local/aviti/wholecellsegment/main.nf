process AVITIWHOLECELLSEGMENT {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    // Reuses the exact container already staged and validated for
    // sopa/segmentationcellpose and CELLPOSEMODEL, which carries Cellpose
    // 4.2.1.1 (cpsam_v2 support) -- no new GPU image needed for the
    // whole-cell path.
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://community.wave.seqera.io/library/python_pip_sopacellpose_cellpose:2bb51160896b005b'
        : 'community.wave.seqera.io/library/python_pip_sopacellpose_cellpose:2bb51160896b005b'}"

    input:
    tuple val(meta), path(nucleus_tif), path(membrane_tif), path(actin_tif)

    output:
    // Named after the source tile, not meta.id: this is exactly the AVITI
    // viewer filename (`<tile>_Cell.tif`) expected under `Well<well>/`, per
    // the Cytocanvas output contract.
    tuple val(meta), path("${meta.tile}_Cell.tif"), emit: cell_mask
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // actin_tif is staged as a literal "NO_FILE" placeholder in 2-channel mode.
    def actin_arg = (actin_tif.name != 'NO_FILE') ? "--actin-tif ${actin_tif}" : ''
    """
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache

    aviti_wholecell_segment.py \\
        --nucleus-tif ${nucleus_tif} \\
        --membrane-tif ${membrane_tif} \\
        ${actin_arg} \\
        --output ${meta.tile}_Cell.tif \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(python3 -c "from importlib.metadata import version; print(version('cellpose'))")
        torch: \$(python3 -c "import torch; print(torch.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tile}_Cell.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(python3 -c "from importlib.metadata import version; print(version('cellpose'))")
        torch: \$(python3 -c "import torch; print(torch.__version__)")
    END_VERSIONS
    """
}
