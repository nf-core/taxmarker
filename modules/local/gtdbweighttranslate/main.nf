process GTDBWEIGHTTRANSLATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(sequences)
    path(bac120_metadata)
    path(ar53_metadata)

    output:
    tuple val(meta), path("*.weights.tsv"), emit: weights
    path "versions.yml",                    emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Nextflow stages an absent optional path as an empty list -- falsy in Groovy.
    def bac120_in = bac120_metadata ? "${bac120_metadata}" : ''
    def ar53_in   = ar53_metadata   ? "${ar53_metadata}"   : ''
    """
    python3 - "${sequences}" "${bac120_in}" "${ar53_in}" "${prefix}.weights.tsv" << 'PYEOF'
import csv
import gzip
import re
import sys
from Bio import SeqIO

sequences_in, bac120_in, ar53_in, weights_out = sys.argv[1:5]

# Same substitution CHECKNAMECONSISTENCY applies to every sequence/taxonomy name --
# duplicated here (not shared) since this runs on the pristine --sequences file,
# before CHECKNAMECONSISTENCY sees it, to recover the header metadata
# (accession, contig length) it hasn't stripped yet.
UNSAFE_CHARS = re.compile(r'[^A-Za-z0-9_.|/-]')

def sanitize(name):
    return UNSAFE_CHARS.sub('_', name)

# GTDB's own ssu_all distribution names each record '<accession>~<contig>', e.g.
# 'RS_GCF_002194975.1~NZ_NHWV01000143.1'. Records without this prefix aren't GTDB
# genome-derived (or are a format this translation doesn't recognise) and are left
# out of the table entirely -- WEIGHTFILTER already treats a missing entry as
# weight 1, so this degrades gracefully on mixed/non-GTDB input.
ACCESSION_RE = re.compile(r'^(RS_|GB_)(GCF|GCA)_[0-9.]+')
CONTIG_LEN_RE = re.compile(r'\\[contig_len=(\\d+)\\]')

RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus']
MISMATCH_PENALTY = {
    'domain': -200, 'phylum': -150, 'class': -40,
    'order': -25, 'family': -10, 'genus': -5,
}
CATEGORY_TERM = {
    'none': 10,                        # isolate
    'derived from single cell': 5,     # SAG
    'derived from metagenome': 0,      # MAG
}
TYPE_MATERIAL_TERM = {
    'na': 0,
    'assembly from type material': 5,
    'assembly from synonym type material': 3,
    'assembly from pathotype material': 3,
    'assembly designated as neotype': 2,
    'assembly designated as reftype': 2,
    'assembly designated as clade exemplar': 2,
}

def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

METADATA_COLUMNS = [
    'ncbi_genome_category', 'gtdb_taxonomy', 'ssu_silva_taxonomy', 'ncbi_type_material_designation',
    'checkm2_completeness', 'checkm2_contamination', 'checkm_completeness', 'checkm_contamination',
]

def sniff_format(path):
    with opener(path) as fh:
        first_line = next((l.strip() for l in fh if l.strip()), '')
    if first_line.startswith('>'):
        return 'fasta'
    if first_line.upper().startswith('CLUSTAL'):
        return 'clustal'
    return 'phylip-relaxed'

with opener(sequences_in) as fh:
    records = list(SeqIO.parse(fh, sniff_format(sequences_in)))
input_accessions = {m.group(0) for m in (ACCESSION_RE.match(r.id) for r in records) if m}

def gtdb_ranks(taxonomy):
    # Strip GTDB's 'x__' rank prefix; only domain..genus (drop species) --
    # the mismatch penalty deliberately stops at genus, see nf-core/taxmarker#15.
    if taxonomy in ('', 'none'):
        return [None] * len(RANKS)
    names = [t.split('__', 1)[-1] for t in taxonomy.split(';')]
    return (names + [None] * len(RANKS))[:len(RANKS)]

def silva_ranks(taxonomy):
    # SILVA taxonomy carries no rank prefix and can be shallower than GTDB's; missing
    # trailing ranks are None, not a mismatch (no evidence either way at that rank).
    if taxonomy in ('', 'none'):
        return [None] * len(RANKS)
    names = taxonomy.split(';')
    return (names + [None] * len(RANKS))[:len(RANKS)]

# --- Build per-rank GTDB-name -> SILVA-synonym sets, from isolate-only co-occurrence
# in this same metadata (see nf-core/taxmarker#15): GTDB and SILVA independently
# revise nomenclature over time, so a real taxonomic synonym (e.g. Desulfobacterota /
# Thermodesulfobacteriota) would otherwise look identical to a genuine contamination
# mismatch. A SILVA name occurring for at least 5% of a GTDB name's isolate rows at
# a given rank is treated as an accepted synonym at that rank (checked empirically
# against real GTDB r226 data: true synonyms sit at ~90-100% co-occurrence, genuine
# errors under ~1%).
SYNONYM_THRESHOLD = 0.05
# Stream the metadata once: synonym counts need every isolate row, but only rows for
# accessions actually in the input are kept, and only the columns used below.
rank_counts = [dict() for _ in RANKS]
metadata = {}
for path in (bac120_in, ar53_in):
    if not path:
        continue
    with opener(path) as fh:
        for row in csv.DictReader(fh, delimiter='\\t'):
            if row['accession'] in input_accessions:
                metadata[row['accession']] = {c: row.get(c) for c in METADATA_COLUMNS}
            if row.get('ncbi_genome_category') != 'none':
                continue
            g_ranks = gtdb_ranks(row.get('gtdb_taxonomy', 'none'))
            s_ranks = silva_ranks(row.get('ssu_silva_taxonomy', 'none'))
            for i, (g_name, s_name) in enumerate(zip(g_ranks, s_ranks)):
                if g_name is None or s_name is None:
                    continue
                counts = rank_counts[i].setdefault(g_name, {})
                counts[s_name] = counts.get(s_name, 0) + 1

synonyms = [dict() for _ in RANKS]
for i, counts_by_gtdb_name in enumerate(rank_counts):
    for g_name, counts in counts_by_gtdb_name.items():
        total = sum(counts.values())
        synonyms[i][g_name] = {s_name for s_name, n in counts.items() if n / total >= SYNONYM_THRESHOLD}

def mismatch_penalty(gtdb_taxonomy, ssu_silva_taxonomy):
    g_ranks = gtdb_ranks(gtdb_taxonomy)
    s_ranks = silva_ranks(ssu_silva_taxonomy)
    for i, rank in enumerate(RANKS):
        g_name, s_name = g_ranks[i], s_ranks[i]
        if g_name is None or s_name is None:
            continue
        if s_name not in synonyms[i].get(g_name, set()):
            return MISMATCH_PENALTY[rank]
    return 0

def quality_score(row):
    completeness, contamination = row.get('checkm2_completeness'), row.get('checkm2_contamination')
    if completeness in (None, '', 'none'):
        completeness, contamination = row.get('checkm_completeness'), row.get('checkm_contamination')
    if completeness in (None, '', 'none') or contamination in (None, '', 'none'):
        return 0.0
    return float(completeness) - 5 * float(contamination)

def contig_len_term(contig_len):
    if contig_len is None:
        return 0
    if contig_len < 2000:
        return 0
    if contig_len < 20000:
        return 1
    if contig_len < 200000:
        return 2
    return 3

def weight_for(row, contig_len):
    return (
        CATEGORY_TERM.get(row.get('ncbi_genome_category'), 0)
        + quality_score(row)
        + TYPE_MATERIAL_TERM.get(row.get('ncbi_type_material_designation'), 0)
        + contig_len_term(contig_len)
        + mismatch_penalty(row.get('gtdb_taxonomy', 'none'), row.get('ssu_silva_taxonomy', 'none'))
    )

with open(weights_out, 'w') as fh:
    for record in records:
        accession_match = ACCESSION_RE.match(record.id)
        if not accession_match:
            continue
        row = metadata.get(accession_match.group(0))
        if row is None:
            continue
        contig_len_match = CONTIG_LEN_RE.search(record.description)
        contig_len = int(contig_len_match.group(1)) if contig_len_match else None
        weight = weight_for(row, contig_len)
        print(f"{sanitize(record.id)}\\t{weight}", file=fh)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.weights.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
