process SATIVAEPANG_MERGEFOLDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sativa-epang:0.10.2--py314hab16a5f_0' :
        'quay.io/biocontainers/sativa-epang:0.10.2--py314hab16a5f_0' }"

    input:
    tuple val(meta), path(taskdirs, stageAs: "shard?/*")

    output:
    tuple val(meta), path("*.l1o_tasks"), emit: taskdir
    tuple val("${task.process}"), val('sativaepang'), eval("sativa-epang --version | cut -d' ' -f2"), topic: versions, emit: versions_sativaepang

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Each shard placed only the folds it was given, so merging every shard's fold
    # contents lands one jplace per fold in any order. Copied, not linked: a link into
    # a shard's work directory dangles when a container mounts only part of work/.
    mkdir "${prefix}.l1o_tasks"
    cp "\$(readlink -f shard1/*/manifest.json)" "${prefix}.l1o_tasks/manifest.json"
    for fold in shard*/*/fold_*; do
        d="${prefix}.l1o_tasks/\$(basename "\$fold")"
        mkdir -p "\$d"
        for f in "\$fold"/*; do
            if [ -e "\$f" ]; then cp -fL "\$f" "\$d/"; fi
        done
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir "${prefix}.l1o_tasks"
    touch "${prefix}.l1o_tasks/manifest.json"
    for fold in shard*/*/fold_*; do
        d="${prefix}.l1o_tasks/\$(basename "\$fold")"
        mkdir -p "\$d"
        touch "\$d/epa_result.jplace"
    done
    """
}
