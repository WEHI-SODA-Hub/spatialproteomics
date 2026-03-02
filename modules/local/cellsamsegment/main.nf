process CELLSAMSEGMENT {
    tag "$meta.id"
    label 'process_multi'
    secret 'DEEPCELL_ACCESS_TOKEN'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/python_pytorch_torchvision_tifffile_pruned:f802da66d91b8999'

    input:
    tuple val(meta), path(tiff), val(nuclear_channel), val(membrane_channels)
    val(compartment)

    output:
    tuple val(meta), path("*.tiff"), emit: segmentation_mask
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_channels = membrane_channels != '' && membrane_channels != [] ? membrane_channels.split(":") : []
    def mem_channel_args = mem_channels.collect { "--membrane-channel '${it}'" }.join(' ')

    """
    cellsam_segment.py \\
        ${tiff} \\
        --output ${prefix}_${compartment}.tiff \\
        --compartment ${compartment} \\
        --nuclear-channel '${nuclear_channel}' \\
        ${mem_channel_args} \\
        --bbox-threshold ${params.cellsam_bbox_threshold} \\
        --block-size ${params.cellsam_block_size} \\
        --overlap ${params.cellsam_overlap} \\
        --iou-depth ${params.cellsam_iou_depth} \\
        --iou-threshold ${params.cellsam_iou_threshold} \\
        ${params.cellsam_use_wsi ? '--use-wsi' : ''} \\
        ${params.cellsam_gauge_cell_size ? '--gauge-cell-size' : ''} \\
        ${params.cellsam_low_contrast_enhancement ? '--low-contrast-enhancement' : ''} \\
        ${params.cellsam_model_path ? "--model-path ${params.cellsam_model_path}" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsam: 0.1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_${compartment}.tiff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsam: 0.1.0
    END_VERSIONS
    """
}
