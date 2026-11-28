process PROFILECOVER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(taxonomy), path(alignment)

    output:
    tuple val(meta), path("*.covfiltered.tax"),   emit: taxonomy
    tuple val(meta), path("*.covfiltered.fasta"), emit: alignment
    tuple val(meta), path("*.cov_excluded.tsv"),  emit: excluded
    path "versions.yml",                          emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ''
    """
    python3 - "${taxonomy}" "${alignment}" \\
        "${prefix}.covfiltered.tax" "${prefix}.covfiltered.fasta" "${prefix}.cov_excluded.tsv" \\
        ${args} << 'PYEOF'
import sys
import argparse
from Bio import SeqIO

# Alignment gap/missing characters -- same set GAPFILTER and RAXTAXFORMAT use.
# Real IUPAC ambiguity codes (e.g. N) count as real content, not gaps: they
# represent an uncertain base call, not an absence of one.
GAP_CHARS = set('-.?')

parser = argparse.ArgumentParser()
parser.add_argument('taxonomy')
parser.add_argument('alignment')
parser.add_argument('out_taxonomy')
parser.add_argument('out_alignment')
parser.add_argument('out_excluded')
parser.add_argument('--min-coverage', type=float, default=0.8,
                     help='Minimum proportion of the HMM profile a sequence must cover '
                          '(non-gap columns in an alignment already masked down to '
                          'match-state columns by HMMER_ESLALIMASK) to be kept; '
                          'sequences below this are too short/incomplete to place '
                          'reliably and are reported separately instead.')
parser.add_argument('--rescue-floor', type=float, default=None,
                     help='For a taxon with no sequence at or above --min-coverage, keep its '
                          'best-covered sequence if that reaches this floor.')
opts = parser.parse_args()

tax_rows = []
with open(opts.taxonomy) as fh:
    for line in fh:
        line = line.rstrip('\\n')
        if not line:
            continue
        name, _, rest = line.partition('\\t')
        tax_rows.append((name, rest))
taxon_of = dict(tax_rows)

records = list(SeqIO.parse(opts.alignment, 'fasta'))

coverage = {}
for record in records:
    seq = str(record.seq)
    non_gap = sum(1 for ch in seq if ch not in GAP_CHARS)
    coverage[record.id] = non_gap / len(seq) if seq else 0.0

kept_names = {name for name, cov in coverage.items() if cov >= opts.min_coverage}
rescued = set()
if opts.rescue_floor is not None:
    covered_taxa = {taxon_of.get(name) for name in kept_names}
    best = {}
    for name, cov in coverage.items():
        taxon = taxon_of.get(name)
        if taxon not in covered_taxa and cov >= opts.rescue_floor and cov > best.get(taxon, ('', -1))[1]:
            best[taxon] = (name, cov)
    rescued = {name for name, _ in best.values()}
    kept_names |= rescued

kept = [record for record in records if record.id in kept_names]
below = [(name, cov) for name, cov in coverage.items() if cov < opts.min_coverage]

with open(opts.out_taxonomy, 'w') as fh:
    for name, rest in tax_rows:
        if name in kept_names:
            print(f"{name}\\t{rest}", file=fh)

SeqIO.write(kept, opts.out_alignment, 'fasta')

with open(opts.out_excluded, 'w') as fh:
    print('seq_name\\tprofile_coverage\\tmin_coverage_threshold\\tstatus', file=fh)
    for name, cov in sorted(below):
        status = 'rescued' if name in rescued else 'excluded'
        print(f"{name}\\t{cov:.4f}\\t{opts.min_coverage}\\t{status}", file=fh)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.covfiltered.tax ${prefix}.covfiltered.fasta ${prefix}.cov_excluded.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
