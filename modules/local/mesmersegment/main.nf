process MESMERSEGMENT {
    tag "$meta.id"
    label 'process_high'
    secret 'DEEPCELL_ACCESS_TOKEN'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/wehi-soda-hub/mesmersegmentation:0.2.0'

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
    def membrane_channel_args = membrane_channels != '' && membrane_channels != [] ?
        membrane_channels.split(":").collect {
            "--membrane-channel \"${it}\""
        }.join(' ') : ''
    """
    mesmer-segment \\
        ${tiff} \\
        --compartment ${compartment} \\
        --nuclear-channel ${nuclear_channel} \\
        ${membrane_channel_args} \\
        ${args} \\
        > "${prefix}_${compartment}.tiff"
    
    # Post-process: transpose mask if dimensions are swapped
    python3 <<'EOF'
import tifffile
import numpy as np

# Read input image to get expected dimensions
input_img = tifffile.imread("${tiff}")
if input_img.ndim == 3:
    expected_shape = input_img.shape[1:]  # (Y, X) from (C, Y, X)
else:
    expected_shape = input_img.shape

# Read mesmer output mask
mask = tifffile.imread("${prefix}_${compartment}.tiff")

# Check if dimensions are swapped (mask is X,Y instead of Y,X)
if mask.shape != expected_shape and mask.shape == expected_shape[::-1]:
    print(f"Transposing mask from {mask.shape} to {expected_shape}")
    mask = mask.T
    tifffile.imwrite("${prefix}_${compartment}.tiff", mask, compression='deflate')
else:
    print(f"Mask dimensions correct: {mask.shape}")
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mesmersegmentation: v0.1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.tiff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mesmersegmentation: v0.1.0
    END_VERSIONS
    """
}
