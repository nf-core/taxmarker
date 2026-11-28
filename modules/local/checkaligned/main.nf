process CHECKALIGNED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.aligned.fasta"),   emit: aligned,   optional: true
    tuple val(meta), path("*.unaligned.fasta"), emit: unaligned, optional: true
    path "versions.yml",                        emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - "${fasta}" "${prefix}.aligned.fasta" "${prefix}.unaligned.fasta" << 'PYEOF'
import shutil
import sys
import Bio
from Bio import SeqIO

fasta_in, aligned_out, unaligned_out = sys.argv[1:4]

lengths = {len(record.seq) for record in SeqIO.parse(fasta_in, 'fasta')}

# A single distinct length across every record means the sequences are already
# aligned (or there's only one record, which is trivially "aligned").
shutil.copy(fasta_in, aligned_out if len(lengths) <= 1 else unaligned_out)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
    print('    biopython: ' + Bio.__version__, file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aligned.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
