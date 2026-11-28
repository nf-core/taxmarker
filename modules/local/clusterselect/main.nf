process CLUSTERSELECT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(taxonomy), path(sequences), path(uc), path(weights)

    output:
    tuple val(meta), path("*.representatives.tax"),   emit: taxonomy
    tuple val(meta), path("*.representatives.fasta"), emit: sequences
    tuple val(meta), path("*.selection.tsv"),         emit: selection
    path "versions.yml",                              emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Nextflow stages an absent optional path(weights) as an empty list -- falsy in
    // Groovy -- distinguishing "no --sequence_weights file given" from a real one.
    def weights_in = weights ? "${weights}" : ''
    """
    python3 - "${taxonomy}" "${sequences}" "${uc}" "${weights_in}" \\
        "${prefix}.representatives.tax" "${prefix}.representatives.fasta" "${prefix}.selection.tsv" << 'PYEOF'
import gzip
import sys
from Bio import SeqIO

taxonomy_in, sequences_in, uc_in, weights_in, out_taxonomy, out_sequences, out_selection = sys.argv[1:8]

# `sequences_in` isn't guaranteed unaligned (see WEIGHTFILTER) -- score on the
# biological (non-gap) length, not raw string length, or an alignment's fixed column
# count would make every candidate the same "length" and the score collapse to weight
# alone.
GAP_CHARS = set('-.?')

# Absent --sequence_weights: every sequence defaults to weight 1 (the multiplicative
# identity) -- see nf-core/taxmarker#15.
weights = {}
if weights_in:
    weights_opener = gzip.open if weights_in.endswith('.gz') else open
    with weights_opener(weights_in, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\\n')
            if not line:
                continue
            seq_id, weight = line.split('\\t')
            weights[seq_id] = float(weight)

taxonomy = {}
with open(taxonomy_in) as fh:
    for line in fh:
        line = line.rstrip('\\n')
        if not line:
            continue
        name, _, lineage = line.partition('\\t')
        taxonomy[name] = lineage

records = {record.id: record for record in SeqIO.parse(sequences_in, 'fasta')}

# vsearch's --uc format: 10 tab-separated columns; S (seed) and H (hit) rows carry
# cluster membership in column 2 (0-based cluster id) and the sequence's own label in
# column 9. C (cluster summary) rows are redundant with the S/H rows and skipped.
cluster_of = {}
opener = gzip.open if uc_in.endswith('.gz') else open
with opener(uc_in, 'rt') as fh:
    for line in fh:
        line = line.rstrip('\\n')
        if not line:
            continue
        fields = line.split('\\t')
        record_type, cluster_id, seq_label = fields[0], fields[1], fields[8]
        if record_type in ('S', 'H'):
            cluster_of[seq_label] = cluster_id

# Score each (cluster, taxon) group's members on length-fraction x weight -- a
# normalized fraction, not raw within-cluster rank, so a cluster's size (driven by
# clustering identity and real population density, not a deliberate choice) can't let
# rank swamp weight's own dynamic range. See nf-core/taxmarker#15.
groups = {}  # (cluster_id, taxon) -> [(seq_id, length, weight, score), ...]
for seq_id, cluster_id in cluster_of.items():
    taxon = taxonomy.get(seq_id)
    if taxon is None:
        continue
    length = sum(1 for ch in str(records[seq_id].seq) if ch not in GAP_CHARS)
    weight = weights.get(seq_id, 1.0)
    groups.setdefault((cluster_id, taxon), []).append([seq_id, length, weight, None])

representative_ids = set()
selection_rows = []
for (cluster_id, taxon), members in groups.items():
    max_length = max(m[1] for m in members)
    for m in members:
        m[3] = (m[1] / max_length) * m[2]
    members.sort(key=lambda m: m[3], reverse=True)
    representative = members[0][0]
    representative_ids.add(representative)
    for seq_id, length, weight, score in members:
        selection_rows.append((cluster_id, taxon, seq_id, length, weight, score, seq_id == representative))

with open(out_taxonomy, 'w') as fh:
    for name in sorted(representative_ids):
        print(f"{name}\\t{taxonomy[name]}", file=fh)

SeqIO.write((records[name] for name in sorted(representative_ids)), out_sequences, 'fasta')

with open(out_selection, 'w') as fh:
    print('cluster_id\\ttaxon\\tseq_name\\tlength\\tweight\\tscore\\tis_representative', file=fh)
    for cluster_id, taxon, seq_id, length, weight, score, is_rep in sorted(selection_rows):
        print(f"{cluster_id}\\t{taxon}\\t{seq_id}\\t{length}\\t{weight:.6f}\\t{score:.6f}\\t{is_rep}", file=fh)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.representatives.tax ${prefix}.representatives.fasta ${prefix}.selection.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
