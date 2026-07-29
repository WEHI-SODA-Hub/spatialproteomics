/*
 * Fetch the Cellpose model weights once per run and emit them as a staged
 * directory, so that segmentation tasks read a file instead of each hitting
 * cellpose.org.
 *
 * Without this, every patch task downloads the ~1.2 GB `cpsam` checkpoint
 * itself: a typical run made ~18 identical requests, which has produced both
 * HTTP 429 rate limiting in CI and truncated downloads locally
 * (`PytorchStreamReader failed reading zip archive`). Staging it once also
 * means `-resume` reuses it across runs, and that every patch task in a run
 * provably reads the same weights.
 *
 * The container is pinned but the weights are not -- cellpose.org serves
 * whatever is current -- so this does not make results reproducible across
 * time. Use `--cellpose_models_dir` with a pre-staged copy for that, which is
 * also what you want on a cluster with no outbound network on compute nodes.
 */
process CELLPOSEMODEL {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.1.11-cellpose'
        : 'docker.io/quentinblampey/sopa:2.1.11-cellpose'}"

    input:
    val _trigger

    output:
    path "cellpose_models", emit: models
    path "versions.yml"   , emit: versions

    script:
    """
    mkdir -p cellpose_models

    # cellpose resolves its model directory at import time, so this must be set
    # before the interpreter starts.
    export CELLPOSE_LOCAL_MODELS_PATH=\$PWD/cellpose_models
    python -c "from cellpose import models; models.CellposeModel(gpu=False)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(python -c "from importlib.metadata import version; print(version('cellpose'))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cellpose_models
    touch cellpose_models/cpsam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(python -c "from importlib.metadata import version; print(version('cellpose'))")
    END_VERSIONS
    """
}
