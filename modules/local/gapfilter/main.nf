process GAPFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(taxonomy), path(alignment)

    output:
    tuple val(meta), path("*.gapfiltered.tax"),   emit: taxonomy
    tuple val(meta), path("*.gapfiltered.fasta"), emit: alignment
    tuple val(meta), path("*.excluded.tsv"),      emit: excluded
    path "versions.yml",                          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ''
    """
    python3 - "${taxonomy}" "${alignment}" \\
        "${prefix}.gapfiltered.tax" "${prefix}.gapfiltered.fasta" "${prefix}.excluded.tsv" \\
        ${args} << 'PYEOF'
import sys
import argparse
from Bio import SeqIO

# Alignment gap/missing characters -- same set RAXTAXFORMAT strips before handing
# sequences to raxtax. Real IUPAC ambiguity codes (e.g. N) count as real content,
# not gaps: they represent an uncertain base call, not an absence of one.
GAP_CHARS = set('-.?')

parser = argparse.ArgumentParser()
parser.add_argument('taxonomy')
parser.add_argument('alignment')
parser.add_argument('out_taxonomy')
parser.add_argument('out_alignment')
parser.add_argument('out_excluded')
parser.add_argument('--min-nongap', type=float, default=0.8,
                     help='Minimum proportion of non-gap columns a sequence must have '
                          'to be kept; sequences below this are too short/incomplete '
                          'to place reliably and are reported separately instead.')
opts = parser.parse_args()

tax_rows = []
with open(opts.taxonomy) as fh:
    for line in fh:
        line = line.rstrip('\\n')
        if not line:
            continue
        name, _, rest = line.partition('\\t')
        tax_rows.append((name, rest))

records = list(SeqIO.parse(opts.alignment, 'fasta'))

kept, excluded = [], []
for record in records:
    seq = str(record.seq)
    non_gap = sum(1 for ch in seq if ch not in GAP_CHARS)
    proportion = non_gap / len(seq) if seq else 0.0
    if proportion >= opts.min_nongap:
        kept.append(record)
    else:
        excluded.append((record.id, proportion))

kept_names = {record.id for record in kept}

with open(opts.out_taxonomy, 'w') as fh:
    for name, rest in tax_rows:
        if name in kept_names:
            print(f"{name}\\t{rest}", file=fh)

SeqIO.write(kept, opts.out_alignment, 'fasta')

with open(opts.out_excluded, 'w') as fh:
    print('seq_name\\tnon_gap_proportion\\tmin_nongap_threshold', file=fh)
    for name, proportion in sorted(excluded):
        print(f"{name}\\t{proportion:.4f}\\t{opts.min_nongap}", file=fh)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gapfiltered.tax ${prefix}.gapfiltered.fasta ${prefix}.excluded.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
