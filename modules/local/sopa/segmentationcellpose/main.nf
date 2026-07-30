/*
 * This module uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */
process SOPA_SEGMENTATIONCELLPOSE {
    // GPU: sopa warns that cellpose >=4 "can be slow without a GPU", and the
    // pipeline now passes --gpu. The reference must stay fully qualified or a
    // configured docker.registry rewrites it and the pull 401s.
    label "process_gpu"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://community.wave.seqera.io/library/python_pip_sopacellpose_cellpose:2bb51160896b005b'
        : 'community.wave.seqera.io/library/python_pip_sopacellpose_cellpose:2bb51160896b005b'}"

    input:
    tuple val(meta), path(zarr), val(index), val(n_patches), val(nuclear_channel), val(membrane_channel)
    path cellpose_models

    output:
    tuple val(meta), path("*.zarr/.sopa_cache/cellpose_boundaries/${index}.parquet"), emit: cellpose_parquet
    path "versions.yml"                                                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def membrane_channel_arg = (membrane_channel && membrane_channel != "[]") ? "--channels \"${membrane_channel}\"" : ""
    """
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache

    # Read the shared weights, but never write into them. ${cellpose_models} is
    # staged as a symlink to a single directory shared by every patch task, so
    # pointing cellpose straight at it would make a download -- if the weights
    # were missing or a new model were requested -- a concurrent write from ~18
    # tasks into one directory. Link the contents into a task-local directory
    # instead: existing weights are still read from the shared copy, and any
    # fetch lands here and harmlessly stays task-local.
    mkdir -p .cellpose_models
    for model in ${cellpose_models}/*; do
        if [ -e "\$model" ]; then
            ln -sf "\$(readlink -f "\$model")" .cellpose_models/
        fi
    done
    export CELLPOSE_LOCAL_MODELS_PATH=\$PWD/.cellpose_models

    sopa segmentation cellpose \\
        ${args} \\
        --patch-index ${index} \\
        ${membrane_channel_arg} \\
        --channels "${nuclear_channel}" \\
        --diameter ${params.cellpose_diameter} \\
        --min-area ${params.cellpose_min_area} \\
        ${zarr}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.zarr/.sopa_cache/cellpose_boundaries
    touch ${prefix}.zarr/.sopa_cache/cellpose_boundaries/${index}.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """
}
