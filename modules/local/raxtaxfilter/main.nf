process RAXTAXFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(raxtax_out), path(taxonomy), path(alignment)

    output:
    tuple val(meta), path("*.filtered.tax"),   emit: taxonomy
    tuple val(meta), path("*.filtered.fasta"), emit: alignment
    tuple val(meta), path("*.mislabels.tsv"),  emit: mislabels
    path "versions.yml",                       emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ''
    """
    python3 - "${raxtax_out}" "${taxonomy}" "${alignment}" \\
        "${prefix}.filtered.tax" "${prefix}.filtered.fasta" "${prefix}.mislabels.tsv" \\
        ${args} << 'PYEOF'
import sys
import argparse
from collections import Counter, defaultdict
from Bio import SeqIO


def load_taxonomy(path):
    tax = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\\n')
            if not line:
                continue
            name, _, lineage = line.partition('\\t')
            tax[name] = [r.strip() for r in lineage.split(';')]
    return tax


def load_raxtax(path):
    # Group every candidate hit by the query's bare sequence name (raxtax reports the
    # full FASTA header, e.g. "SeqA;tax=...;", as the query identifier for a self-vs-self
    # run since query and database are the same file).
    hits = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\\n')
            if not line:
                continue
            query, lineage, confidences, _local_signal, _global_signal = line.split('\\t')
            bare_name = query.split(';', 1)[0]
            predicted = lineage.split(',')
            conf_values = [float(c) for c in confidences.split(',')]
            hits[bare_name].append((predicted, conf_values))
    return hits


parser = argparse.ArgumentParser()
parser.add_argument('raxtax_out')
parser.add_argument('taxonomy')
parser.add_argument('alignment')
parser.add_argument('out_taxonomy')
parser.add_argument('out_alignment')
parser.add_argument('out_mislabels')
parser.add_argument('--filter-rank', type=int, default=1,
                     help='Leaf-counted rank (1 = most specific) at which the best raxtax '
                          'hit must agree with the declared taxonomy; disagreement there '
                          'flags the sequence as a likely mislabel.')
opts = parser.parse_args()

tax = load_taxonomy(opts.taxonomy)
hits = load_raxtax(opts.raxtax_out)

# --skip-exact-matches hides a query from its own classification, so a taxon with no
# other member can never be predicted and would always look mislabelled.
members = Counter(tuple(lineage[:i + 1]) for lineage in tax.values() for i in range(len(lineage)))

flagged = {}
for name, candidates in hits.items():
    declared = tax.get(name)
    if declared is None:
        continue

    # Best hit = the candidate with the highest confidence at its deepest (most
    # specific) reported rank -- the single most likely classification for this
    # query, as opposed to SATIVASCORE's weighted aggregate across all candidates.
    best_predicted, best_conf = max(candidates, key=lambda c: c[1][-1])

    n_ranks = min(len(declared), len(best_predicted))
    root_index = n_ranks - opts.filter_rank  # 0-indexed, root-counted position to check
    if root_index < 0:
        continue  # requested rank goes deeper than this lineage; nothing to check
    if members[tuple(declared[:root_index + 1])] < 2:
        continue

    if declared[root_index] != best_predicted[root_index]:
        confidence = best_conf[root_index] if root_index < len(best_conf) else best_conf[-1]
        flagged[name] = {
            'original': declared,
            'predicted': best_predicted,
            'confidence': confidence,
            'mismatch_rank': root_index + 1,
        }

kept_names = set(tax) - set(flagged)

with open(opts.out_taxonomy, 'w') as fh:
    for name in sorted(kept_names):
        print(f"{name}\\t{';'.join(tax[name])}", file=fh)

with open(opts.out_alignment, 'w') as fh:
    for record in SeqIO.parse(opts.alignment, 'fasta'):
        if record.id in kept_names:
            print(f">{record.id}", file=fh)
            print(str(record.seq), file=fh)

with open(opts.out_mislabels, 'w') as fh:
    print('seq_name\\toriginal_label\\tpredicted_label\\tlwr\\tmismatch_rank\\tmethod', file=fh)
    for name in sorted(flagged):
        info = flagged[name]
        print(f"{name}\\t{';'.join(info['original'])}\\t{';'.join(info['predicted'])}\\t"
              f"{info['confidence']:.6f}\\t{info['mismatch_rank']}\\traxtax", file=fh)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.tax ${prefix}.filtered.fasta ${prefix}.mislabels.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
