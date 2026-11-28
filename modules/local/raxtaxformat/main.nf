process RAXTAXFORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(alignment), path(taxonomy)

    output:
    tuple val(meta), path("*.raxtax.fasta"), emit: fasta
    path "versions.yml",                     emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - "${alignment}" "${taxonomy}" "${prefix}.raxtax.fasta" << 'PYEOF'
import sys
from Bio import SeqIO

alignment, taxonomy, out_fasta = sys.argv[1:4]

# Alignment gap/missing characters to strip so k-mer classification isn't
# disrupted by columns that carry no sequence information. Real IUPAC
# ambiguity codes (e.g. N) are kept: they represent an uncertain base call,
# not an alignment gap. U/u (RNA uracil) is normalised to T/t: raxtax's
# parser only understands DNA and panics on 'U' (observed on this pipeline's
# own RNA test data), but U and T represent the same base for classification.
SEQ_TRANSLATION = str.maketrans({c: None for c in '-.?'} | {'U': 'T', 'u': 't'})

tax = {}
with open(taxonomy) as fh:
    for line in fh:
        line = line.rstrip('\\n')
        if not line:
            continue
        name, _, lineage = line.partition('\\t')
        tax[name] = lineage

with open(out_fasta, 'w') as out_fh:
    for record in SeqIO.parse(alignment, 'fasta'):
        lineage = tax.get(record.id, '')
        raxtax_tax = ','.join(part.strip() for part in lineage.split(';'))
        seq = str(record.seq).translate(SEQ_TRANSLATION)
        print(f'>{record.id};tax={raxtax_tax};', file=out_fh)
        print(seq, file=out_fh)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raxtax.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
