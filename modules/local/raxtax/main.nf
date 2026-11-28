process RAXTAX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raxtax:1.6.0--h4349ce8_0' :
        'quay.io/biocontainers/raxtax:1.6.0--h4349ce8_0' }"

    input:
    tuple val(meta), path(query)
    path(database, stageAs: 'database/*')

    output:
    tuple val(meta), path("${prefix}/raxtax.out"),           emit: out
    tuple val(meta), path("${prefix}/raxtax.log"),           emit: log
    tuple val(meta), path("${prefix}/raxtax.tsv"), optional: true, emit: tsv
    path "versions.yml",                                     emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    raxtax \\
        -d ${database} \\
        -i ${query} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxtax: \$(raxtax --version | sed 's/^raxtax //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/raxtax.out
    touch ${prefix}/raxtax.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxtax: \$(raxtax --version | sed 's/^raxtax //')
    END_VERSIONS
    """
}
