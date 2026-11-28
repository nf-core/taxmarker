process EXPORTREFERENCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(sequences), path(taxonomy), path(selection), path(verified, stageAs: 'verified/*'), path(mislabels)
    val(counts)

    output:
    tuple val(meta), path("*.addSpecies.fna.gz"),     emit: addspecies
    tuple val(meta), path("*.assignTaxonomy.fna.gz"), emit: assigntaxonomy
    tuple val(meta), path("*.general.fna.gz"),        emit: general
    tuple val(meta), path("*.names.txt"),             emit: names
    path "versions.yml",                              emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Absent optional paths are staged as empty lists, falsy in Groovy.
    def selection_in = selection ? "${selection}" : ''
    def mislabels_in = mislabels ? "${mislabels}" : ''
    """
    python3 - "${sequences}" "${taxonomy}" "${selection_in}" "${verified}" "${mislabels_in}" \\
        "${prefix}" ${counts.join(' ')} << 'PYEOF'
import csv
import gzip
import re
import sys
from Bio import SeqIO

sequences_in, taxonomy_in, selection_in, verified_in, mislabels_in, prefix = sys.argv[1:7]
counts = [None if n == 'all' else int(n) for n in sys.argv[7:]]

# Must match CHECKNAMECONSISTENCY, which renamed every sequence the other inputs name.
UNSAFE_CHARS = re.compile(r'[^A-Za-z0-9_.|/-]')
GAP_CHARS = '-.?'

def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

def sniff_format(path):
    with opener(path) as fh:
        first_line = next((l.strip() for l in fh if l.strip()), '')
    if first_line.startswith('>'):
        return 'fasta'
    if first_line.upper().startswith('CLUSTAL'):
        return 'clustal'
    return 'phylip-relaxed'

def read_names(path):
    with open(path) as fh:
        return {line.split('\\t', 1)[0] for line in fh if line.strip()}

records = {}
with opener(sequences_in) as fh:
    for record in SeqIO.parse(fh, sniff_format(sequences_in)):
        seq = str(record.seq)
        for ch in GAP_CHARS:
            seq = seq.replace(ch, '')
        records[UNSAFE_CHARS.sub('_', record.id)] = (record.description, seq)

verified = read_names(verified_in)
flagged = set()
if mislabels_in:
    with open(mislabels_in) as fh:
        flagged = {row['seq_name'] for row in csv.DictReader(fh, delimiter='\\t')}

# (name, taxon, length, weight, representative)
rows = []
if selection_in:
    with open(selection_in) as fh:
        selection = list(csv.DictReader(fh, delimiter='\\t'))
    rep_of = {(r['cluster_id'], r['taxon']): r['seq_name'] for r in selection if r['is_representative'] == 'True'}
    for r in selection:
        rows.append((r['seq_name'], r['taxon'], int(r['length']), float(r['weight']), rep_of[(r['cluster_id'], r['taxon'])]))
else:
    with open(taxonomy_in) as fh:
        for line in fh:
            if line.strip():
                name, taxon = line.rstrip('\\n').split('\\t', 1)
                rows.append((name, taxon, len(records[name][1]), 1.0, name))

# Rounded weights can lift a member above the representative CLUSTERSELECT chose on the
# unrounded values, so a member never scores above its representative.
score = {name: round(weight * length, 3) for name, _, length, weight, _ in rows}
by_taxon = {}
for name, taxon, length, weight, rep in rows:
    if rep in verified and rep not in flagged:
        by_taxon.setdefault(taxon, []).append((-min(score[name], score[rep]), name != rep, name))
for candidates in by_taxon.values():
    candidates.sort()

def labels(taxon):
    return [re.sub(r'^[a-z]__', '', rank) for rank in taxon.split(';')]

# Without rank prefixes the last two ranks are taken as genus and species, which is only
# right for a lineage as deep as the deepest one.
max_depth = max((len(taxon.split(';')) for taxon in by_taxon), default=0)

def species_name(taxon):
    ranks = taxon.split(';')
    if any(re.match(r'^[a-z]__', r) for r in ranks):
        by_rank = {r[0]: r[3:] for r in ranks if re.match(r'^[a-z]__', r)}
        genus, species = by_rank.get('g', ''), by_rank.get('s', '')
    elif len(ranks) >= max(2, max_depth):
        genus, species = ranks[-2], ranks[-1]
    else:
        return None
    if not genus.strip() or not species.strip():
        return None
    return species if species.startswith(genus + ' ') else f"{genus} {species}"

for n in counts:
    label = 'all' if n is None else f"n{n}"
    with gzip.open(f"{prefix}.{label}.general.fna.gz", 'wt') as general, \
         gzip.open(f"{prefix}.{label}.addSpecies.fna.gz", 'wt') as add_species, \
         gzip.open(f"{prefix}.{label}.assignTaxonomy.fna.gz", 'wt') as assign_taxonomy, \
         open(f"{prefix}.{label}.names.txt", 'w') as names:
        for taxon in sorted(by_taxon):
            for *_, name in by_taxon[taxon][:n]:
                description, seq = records[name]
                seq_id = description.split()[0]
                print(f">{description}\\n{seq}", file=general)
                if species_name(taxon):
                    print(f">{seq_id} {species_name(taxon)}\\n{seq}", file=add_species)
                # Trailing ';': DADA2 appends the padding for a shorter lineage straight onto the header.
                print(f">{';'.join(labels(taxon))};\\n{seq}", file=assign_taxonomy)
                print(name, file=names)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def files = counts.collectMany { n -> ["addSpecies", "assignTaxonomy", "general"].collect { f -> "${prefix}.${n == 'all' ? 'all' : 'n' + n}.${f}.fna.gz" } }.join(' ')
    def names = counts.collect { n -> "${prefix}.${n == 'all' ? 'all' : 'n' + n}.names.txt" }.join(' ')
    """
    for f in ${files}; do echo -n | gzip > \$f; done
    touch ${names}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
