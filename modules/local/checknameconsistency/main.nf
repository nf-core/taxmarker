
process CHECKNAMECONSISTENCY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(taxonomy), path(alignment)

    output:
    tuple val(meta), path("*.checked.tax"), path("*.checked.${alignment.extension}"), emit: checked
    path "versions.yml",                                                             emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Output names must differ from the input taxonomy/alignment filenames: Nextflow
    // stages inputs as symlinks, and Python's open(path, 'w') follows a symlink rather
    // than replacing it, so reusing the input's exact name would silently overwrite the
    // original source file through the symlink.
    """
    python3 - "${taxonomy}" "${alignment}" "${prefix}.checked.tax" "${prefix}.checked.${alignment.extension}" << 'PYEOF'
import re
import sys
import Bio
from Bio import SeqIO

taxonomy_in, alignment_in, taxonomy_out, alignment_out = sys.argv[1:5]

# IQTREE's own tree-writer silently replaces almost every non-alphanumeric
# character with '_' when writing leaf names (empirically confirmed: parens,
# '~', '#', and everything else tried except '|' and '/' get mangled this
# way). Left alone here, that would desync alignment/taxonomy names from the
# tree's leaf names downstream, so match IQTREE's own behaviour up front
# rather than chasing individual problem characters as they turn up.
UNSAFE_CHARS = re.compile(r'[^A-Za-z0-9_.|/-]')

def sanitize(name):
    return UNSAFE_CHARS.sub('_', name)

def find_duplicates(names):
    seen, dups = set(), set()
    for name in names:
        if name in seen:
            dups.add(name)
        seen.add(name)
    return dups

# --- Taxonomy: seq_name<TAB>rank1;rank2;... ---
tax_rows = []
with open(taxonomy_in) as fh:
    for line in fh:
        line = line.rstrip('\\n')
        if not line:
            continue
        name, _, rest = line.partition('\\t')
        tax_rows.append((sanitize(name), rest))
tax_names = [name for name, _ in tax_rows]

# --- Alignment: format sniffed from content (FASTA, Clustal or PHYLIP; this
# runs before the pipeline's own EMBOSS_SEQRET normalisation to PHYLIP) ---
with open(alignment_in) as fh:
    first_line = next((l.strip() for l in fh if l.strip()), '')
if first_line.startswith('>'):
    alignment_format = 'fasta'
elif first_line.upper().startswith('CLUSTAL'):
    alignment_format = 'clustal'
else:
    alignment_format = 'phylip-relaxed'

records = list(SeqIO.parse(alignment_in, alignment_format))
for record in records:
    record.id = sanitize(record.id)
    record.name = record.id
    record.description = record.id
aln_names = [record.id for record in records]

# --- Validate: no post-rewrite collisions, same name set in both files ---
problems = []

tax_dups = find_duplicates(tax_names)
if tax_dups:
    problems.append('Duplicate names in taxonomy after rewriting problematic characters: ' + ', '.join(sorted(tax_dups)))

aln_dups = find_duplicates(aln_names)
if aln_dups:
    problems.append('Duplicate names in alignment after rewriting problematic characters: ' + ', '.join(sorted(aln_dups)))

only_in_alignment = sorted(set(aln_names) - set(tax_names))
only_in_taxonomy  = sorted(set(tax_names) - set(aln_names))
if only_in_alignment:
    problems.append('Names in alignment but not in taxonomy: ' + ', '.join(only_in_alignment))
if only_in_taxonomy:
    problems.append('Names in taxonomy but not in alignment: ' + ', '.join(only_in_taxonomy))

if problems:
    sys.exit('\\n'.join(problems))

# --- Write out the (possibly rewritten) files ---
with open(taxonomy_out, 'w') as fh:
    for name, rest in tax_rows:
        print(f"{name}\\t{rest}", file=fh)

SeqIO.write(records, alignment_out, alignment_format)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
    print('    biopython: ' + Bio.__version__, file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.checked.tax ${prefix}.checked.${alignment.extension}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
