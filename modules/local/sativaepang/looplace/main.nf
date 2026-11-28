process SATIVAEPANG_LOOPLACE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sativa-epang:0.10.2--py314hab16a5f_0' :
        'quay.io/biocontainers/sativa-epang:0.10.2--py314hab16a5f_0' }"

    input:
    tuple val(meta), path(taskdir, stageAs: "input")

    output:
    tuple val(meta), path("*.l1o_tasks"), emit: taskdir
    tuple val("${task.process}"), val('sativaepang'), eval("sativa-epang --version | cut -d' ' -f2"), topic: versions, emit: versions_sativaepang

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # loo-place writes into every fold, so it works in a private directory, not the
    # read-only staged input. The fold inputs are linked, not copied (near-full alignments
    # per fold), and the links removed afterwards: looscore reads only manifest.json and
    # the placements, and a link into lootasks' work directory dangles when a container
    # mounts only part of work/.
    mkdir "${prefix}.l1o_tasks"
    cp "\$(readlink -f input/manifest.json)" "${prefix}.l1o_tasks/manifest.json"
    for fold in input/fold_*; do
        d="${prefix}.l1o_tasks/\$(basename "\$fold")"
        mkdir "\$d"
        ln -s "\$(readlink -f "\$fold")"/* "\$d/"
    done

    sativa-epang \\
        -stage loo-place \\
        -taskdir "${prefix}.l1o_tasks" \\
        -T ${task.cpus} \\
        ${args}

    find "${prefix}.l1o_tasks" -type l -exec rm -f {} +
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir "${prefix}.l1o_tasks"
    cp "\$(readlink -f input/manifest.json)" "${prefix}.l1o_tasks/manifest.json"
    for fold in input/fold_*; do
        d="${prefix}.l1o_tasks/\$(basename "\$fold")"
        mkdir "\$d"
        touch "\$d/epa_result.jplace"
    done
    """
}
