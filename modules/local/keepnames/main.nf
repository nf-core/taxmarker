process KEEPNAMES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(taxonomy), path(alignment), path(names)

    output:
    tuple val(meta), path("*.reference.tax"),   emit: taxonomy
    tuple val(meta), path("*.reference.fasta"), emit: alignment
    path "versions.yml",                        emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - "${taxonomy}" "${alignment}" "${names}" "${prefix}.reference.tax" "${prefix}.reference.fasta" << 'PYEOF'
import sys

taxonomy_in, alignment_in, names_in, taxonomy_out, alignment_out = sys.argv[1:6]

with open(names_in) as fh:
    wanted = {line.strip() for line in fh if line.strip()}

with open(taxonomy_in) as fh, open(taxonomy_out, 'w') as out:
    for line in fh:
        if line.strip() and line.split('\\t', 1)[0] in wanted:
            out.write(line)

keep = False
seen = set()
with open(alignment_in) as fh, open(alignment_out, 'w') as out:
    for line in fh:
        if line.startswith('>'):
            name = line[1:].split()[0]
            keep = name in wanted
            if keep:
                seen.add(name)
        if keep:
            out.write(line)

if seen != wanted:
    sys.exit('Not in the alignment: ' + ', '.join(sorted(wanted - seen)))

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reference.tax ${prefix}.reference.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
