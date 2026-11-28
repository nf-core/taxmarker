process WEIGHTFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(taxonomy), path(sequences), path(weights), val(min_weight)

    output:
    tuple val(meta), path("*.weightfiltered.tax"),         emit: taxonomy
    tuple val(meta), path("*.weightfiltered.fasta"),       emit: sequences
    tuple val(meta), path("*.weightfiltered.ungapped.fasta"), emit: ungapped_sequences
    tuple val(meta), path("*.excluded.tsv"),               emit: excluded
    path "versions.yml",                                   emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Nextflow stages an absent optional path(weights) as an empty list -- falsy in
    // Groovy -- distinguishing "no --sequence_weights file given" from a real one.
    def weights_in = weights ? "${weights}" : ''
    """
    python3 - "${taxonomy}" "${sequences}" "${weights_in}" "${min_weight}" \\
        "${prefix}.weightfiltered.tax" "${prefix}.weightfiltered.fasta" "${prefix}.weightfiltered.ungapped.fasta" "${prefix}.excluded.tsv" << 'PYEOF'
import gzip
import sys
from Bio import SeqIO
from Bio.Seq import Seq

taxonomy_in, sequences_in, weights_in, min_weight, out_taxonomy, out_sequences, out_ungapped, out_excluded = sys.argv[1:9]
min_weight = float(min_weight)

def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

# The caller doesn't guarantee unaligned input at this point (ENSURE_ALIGNED, which
# detects aligned-vs-not, runs later) -- VSEARCH_CLUSTER needs plain sequences, same
# gap/missing character set RAXTAXFORMAT strips for raxtax. Emitted as a SEPARATE
# ungapped copy alongside the untouched original: degapping every kept record in place
# would make already-aligned input look unaligned to ENSURE_ALIGNED's own length-based
# detection once representatives are picked back out downstream.
GAP_CHARS = str.maketrans({c: None for c in '-.?'})

# Absent --sequence_weights: every sequence defaults to weight 1 (the multiplicative
# identity), not an imposed scale -- see nf-core/taxmarker#15.
weights = {}
if weights_in:
    with opener(weights_in) as fh:
        for line in fh:
            line = line.rstrip('\\n')
            if not line:
                continue
            seq_id, weight = line.split('\\t')
            weights[seq_id] = float(weight)

tax_rows = []
with open(taxonomy_in) as fh:
    for line in fh:
        line = line.rstrip('\\n')
        if not line:
            continue
        name, _, rest = line.partition('\\t')
        tax_rows.append((name, rest))

records = list(SeqIO.parse(sequences_in, 'fasta'))

kept, excluded = [], []
for record in records:
    weight = weights.get(record.id, 1.0)
    if weight >= min_weight:
        kept.append(record)
    else:
        excluded.append((record.id, weight))

kept_names = {record.id for record in kept}

with open(out_taxonomy, 'w') as fh:
    for name, rest in tax_rows:
        if name in kept_names:
            print(f"{name}\\t{rest}", file=fh)

SeqIO.write(kept, out_sequences, 'fasta')

def ungapped(record):
    # Bio.Seq.Seq.translate() means codon->protein translation, not str.translate() --
    # go through str explicitly to just strip gap characters. Copy the record rather
    # than mutating it, since `kept` (with gaps intact) was already written above.
    ungapped_record = record[:]
    ungapped_record.seq = Seq(str(record.seq).translate(GAP_CHARS))
    return ungapped_record

SeqIO.write((ungapped(record) for record in kept), out_ungapped, 'fasta')

with open(out_excluded, 'w') as fh:
    print('seq_name\\tweight\\tmin_weight_threshold', file=fh)
    for name, weight in sorted(excluded):
        print(f"{name}\\t{weight:.6f}\\t{min_weight}", file=fh)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.weightfiltered.tax ${prefix}.weightfiltered.fasta ${prefix}.weightfiltered.ungapped.fasta ${prefix}.excluded.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
