process CREATEREPORT {
    tag "$meta.id"
    label 'process_single'
    label 'process_medium_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'docker://ghcr.io/wehi-soda-hub/spatialvis:latest' :
                'ghcr.io/wehi-soda-hub/spatialvis:latest' }"
    input:
    tuple val(meta),
        path(expression_file),
        path(hierarchy_file),
        val(sample_id),
        val(marker_col),
        val(markers),
        val(are_markers_split),
        val(cell_types),
        val(parent_types),
        val(metadata_cols),
        val(plot_heatmaps),
        val(plot_props),
        val(plot_umap),
        val(plot_clusters),
        val(plot_spatial),
        val(save_rdata)
    val(template_file)

    output:
    tuple val(meta), path("*/*.html"), path("*/*_files/*"), emit: report
    tuple val(meta), path("*/*.rds"), emit: rds, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    if [[ ! -f report_template.qmd ]]; then
        ln -s ${template_file} report_template.qmd
    fi
    quarto render report_template.qmd \\
        --to html \\
        --no-cache \\
        --output ${prefix}.html \\
        ${args} \\
        -P hierarchy_file:${hierarchy_file} \\
        -P expression_file:${expression_file} \\
        -P sample_id:${sample_id} \\
        -P marker_col:${marker_col} \\
        -P markers:${markers} \\
        -P are_markers_split:${are_markers_split} \\
        -P cell_types:${cell_types} \\
        -P parent_types:${parent_types} \\
        -P metadata_cols:${metadata_cols} \\
        -P plot_heatmaps:${plot_heatmaps} \\
        -P plot_props:${plot_props} \\
        -P plot_umap:${plot_umap} \\
        -P plot_clusters:${plot_clusters} \\
        -P plot_spatial:${plot_spatial} \\
        -P save_rdata:${save_rdata} \\
        -P sample_name:${meta.id}

    mkdir -p ${prefix}
    mv ${prefix}.html ${prefix}
    if [[ -f ${prefix}.rds ]]; then
        mv ${prefix}.rds ${prefix}
    fi
    cp -r report_template_files ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        spatialVis: \$(Rscript -e "cat(as.character(packageVersion('spatialVis')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/report_template_files
    touch ${prefix}/${prefix}.html
    touch ${prefix}/${prefix}.rds
    touch ${prefix}/report_template_files/foo.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        spatialVis: \$(Rscript -e "cat(as.character(packageVersion('spatialVis')))")
    END_VERSIONS
    """
}
