process TREEGRAFT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(plan), path(results), path(constraint)

    output:
    tuple val(meta), path("*.grafted.nwk"), emit: tree
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    treedivide.py graft \\
        --plan ${plan} \\
        --constraint ${constraint} \\
        --out ${prefix}.grafted.nwk \\
        ${args} \\
        ${results}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '();' > ${prefix}.grafted.nwk
    """
}
