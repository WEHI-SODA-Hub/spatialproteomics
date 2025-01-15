process CREATEREPORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'docker://ghcr.io/wehi-soda-hub/spatialvis:latest' :
                'ghcr.io/wehi-soda-hub/spatialvis:latest' }"
    input:
    tuple val(meta), path(expression_file), path(hierarchy_file), val(markers), val(cell_types), val(parent_types), val(downstream_analyses)
    path(template_file)

    output:
    tuple val(meta), path("*.html"), path("*_files/*"), emit: report
    path("*.rds"), emit: rds
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    quarto render ${template_file} \\
        --to html \\
        --no-cache \\
        --output ${prefix}.html \\
        -P hierarchies_file:${hierarchy_file} \\
        -P expression_file:${expression_file} \\
        -P sample_name:${meta.id}

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
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        spatialVis: \$(Rscript -e "cat(as.character(packageVersion('spatialVis')))")
    END_VERSIONS
    """
}
