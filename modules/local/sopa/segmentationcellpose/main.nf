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
    def has_membrane = (membrane_channel && membrane_channel != "[]")
    def membrane_channel_arg = has_membrane ? "--channels \"${membrane_channel}\"" : ""

    // Cellpose's normalize dict is single-sourced here (rather than in
    // modules.config) because global normalisation needs the per-sample bounds
    // that only exist at run time. sopa forwards --method-kwargs straight into
    // model.eval, and passing it twice would conflict, so exactly one
    // method_kwargs.json is written below and always passed.
    //
    // Precedence: global normalisation > tile-local percentile > cellpose's
    // stock whole-patch percentile ({} == no override).
    def tile_norm = (params.containsKey('cellpose_tile_norm_blocksize') && params.cellpose_tile_norm_blocksize) ? params.cellpose_tile_norm_blocksize : 0
    def global_norm = params.containsKey('cellpose_normalize_global') && params.cellpose_normalize_global
    def percentiles = (params.containsKey('cellpose_normalize_percentiles') && params.cellpose_normalize_percentiles) ? params.cellpose_normalize_percentiles : '1,99'
    // Channels in the SAME order sopa builds the image below (membrane then
    // nucleus), so the per-channel bounds line up with the array cellpose sees.
    def norm_channel_args = (has_membrane ? "--channel \"${membrane_channel}\" " : "") + "--channel \"${nuclear_channel}\""
    def norm_setup = ''
    if (global_norm) {
        norm_setup = "compute_cellpose_norm.py ${zarr} ${norm_channel_args} --percentiles ${percentiles} --as-method-kwargs -o method_kwargs.json"
    } else if (tile_norm) {
        norm_setup = "printf '%s' '{\"normalize\": {\"tile_norm_blocksize\": ${tile_norm}, \"percentile\": [1.0, 99.0]}}' > method_kwargs.json"
    } else {
        norm_setup = "printf '%s' '{}' > method_kwargs.json"
    }
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

    # Write exactly one cellpose normalize dict (see script block above).
    ${norm_setup}

    sopa segmentation cellpose \\
        ${args} \\
        --patch-index ${index} \\
        ${membrane_channel_arg} \\
        --channels "${nuclear_channel}" \\
        --diameter ${params.cellpose_diameter} \\
        --min-area ${params.cellpose_min_area} \\
        --method-kwargs "\$(cat method_kwargs.json)" \\
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
