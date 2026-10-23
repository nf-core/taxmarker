process TREEDIVIDE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(alignment), path(constraint), path(partition_tree)

    output:
    tuple val(meta), path("parts/*.fasta"),         emit: fasta       , optional: true
    tuple val(meta), path("parts/*.nwk"),           emit: constraints , optional: true
    tuple val(meta), path("parts/resolved/*.nwk"),  emit: resolved    , optional: true
    tuple val(meta), path("*.plan.json"),    emit: plan
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def ptree_in = partition_tree ? "--partition-tree ${partition_tree}" : ''
    """
    treedivide.py split \\
        --alignment ${alignment} \\
        --constraint ${constraint} \\
        ${ptree_in} \\
        --outdir parts \\
        --plan ${prefix}.plan.json \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir parts
    touch parts/job_000.fasta parts/job_000.nwk
    echo '{"jobs": [{"id": "job_000", "n": 0, "reps": [], "outgroups": [], "root": true, "resolved": false}]}' > ${prefix}.plan.json
    """
}
