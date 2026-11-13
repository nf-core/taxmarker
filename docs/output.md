# nf-core/taxmarker: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

The final products of a run are in two places: the table of putative mislabels in [`mislabels/`](#mislabels), and everything meant for use afterwards in [`export/`](#export).
The other directories hold intermediate results, one directory per pipeline step, in the order in which the steps run.
They are mostly useful to follow what happened to a sequence, and most users can skip them.

A directory is only created when its step runs: some steps depend on the input (for example, alignment only runs for unaligned sequences) and some are switched on or off by parameters.
Each section below says when its directory appears.

Most files are named after `user-alignment`, the name the pipeline gives to the input set.
Taxonomy files (`*.tax`) are tab-separated, without a header, with two columns: the sequence name and its lineage as a `;`-separated list of ranks, for example `Bacteria;Firmicutes;Clostridia`.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- Final results
  - [Mislabels](#mislabels) - Putative mislabels found by sativa-epang and raxtax
  - [Export](#export) - Final products for use afterwards: a copy of the mislabels table, reference sequences and a reference phylogeny
- Input preparation
  - [Sequence selection](#sequence-selection) - Keep only the sequences whose header matches `--seqgrep`
  - [Taxonomy resolution](#taxonomy-resolution) - Read the taxonomy from `--taxonomy` or from the sequence headers
  - [Name check](#name-check) - Check that sequences and taxonomy name the same sequences, and make the names safe
  - [Conversion to FASTA](#conversion-to-fasta) - Convert the sequences to FASTA with EMBOSS seqret
- Weighted clustering
  - [GTDB weights](#gtdb-weights) - Sequence weights computed from GTDB genome metadata
  - [Weight filter](#weight-filter) - Drop sequences below `--min_weight` and make an ungapped copy for VSEARCH
  - [Clustering](#clustering) - Cluster identical or similar sequences with VSEARCH
  - [Cluster selection](#cluster-selection) - Keep one representative per cluster and taxon
- raxtax prefilter
  - [raxtax input](#raxtax-input) - Sequences in the format raxtax reads
  - [raxtax classification](#raxtax-classification) - Classify every sequence against all the others
  - [raxtax filter](#raxtax-filter) - Flag sequences that raxtax confidently places in another taxon
- Alignment
  - [Alignment check](#alignment-check) - Find out whether the sequences are already aligned
  - [Alignment](#alignment) - Align unaligned sequences with hmmalign or MAFFT
  - [HMM profile](#hmm-profile) - The profile picked from `--hmm` by `--hmm_name`
- Short-sequence filters
  - [Gap filter](#gap-filter) - Drop sequences with too many gaps (aligned input)
  - [Profile coverage filter](#profile-coverage-filter) - Drop sequences that cover too little of the HMM profile (hmmalign input)
- Reference tree and placement
  - [Taxonomy tree](#taxonomy-tree) - The taxonomy as a tree, used as a constraint for the tree search
  - [Reference tree](#reference-tree) - The reference phylogeny built by RAxML-NG
  - [Leave-one-out placement](#leave-one-out-placement) - sativa-epang's reference, placements and scores
- Reporting
  - [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
  - [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Final results

### Mislabels

<details markdown="1">
<summary>Output files</summary>

- `mislabels/`
  - `*.mislabels.tsv`: putative mislabels, one row per flagged sequence, from sativa-epang and from the raxtax prefilter.
- `export/`
  - `*.mislabels.tsv`: a copy of the same table, next to the other final products.

</details>

The table joins the rows of the two methods: the raxtax rows come from [`raxtaxfilter/`](#raxtax-filter) and the sativa rows from [`sativaepang/`](#leave-one-out-placement).
A sequence flagged by raxtax does not go on to the placement step, so it is never in both.

The table is tab-separated, with a header line and these columns:

| Column                     | Description                                                                                                                                |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `seq_name`                 | Sequence name                                                                                                                              |
| `original_label`           | Declared taxonomy, as a `;`-separated lineage                                                                                              |
| `predicted_label`          | Taxonomy proposed by the method                                                                                                            |
| `lwr`                      | Confidence at the mismatching rank: sativa-epang's `Confidence`, or the raxtax confidence                                                  |
| `original_label_exclusion` | sativa-epang's `OriginalLabelExclusion`: the placement weight outside the declared taxon at the mislabelled rank. `NA` for raxtax rows     |
| `mismatch_rank`            | Rank where the lineages disagree, counted from the root (1 = first rank): the first differing rank for sativa, the checked rank for raxtax |
| `method`                   | `sativa` or `raxtax`                                                                                                                       |

A sativa-epang row with `lwr` 0.51 and `original_label_exclusion` near 1 is a stem placement: the sequence was placed on the branch between the proposed taxon and its parent.
Sativa splits the weight of such a placement 0.49/0.51 between the two (a hardcoded value), so the 0.51 does not come from the data.
The sequence is confidently outside its declared taxon, but it is not clear which neighbouring taxon it belongs to.

### Export

<details markdown="1">
<summary>Output files</summary>

- `export/`
  - `*.mislabels.tsv`: a copy of the [mislabels table](#mislabels). Written whenever raxtax or the placement step runs.
  - With `--export_n_per_species`, four files for each number N (`--run_reference` also writes them, for N = 1):
    - `*.n<N>.addSpecies.fna.gz`: reference for DADA2's `addSpecies()`, headers `>ID Genus species`.
    - `*.n<N>.assignTaxonomy.fna.gz`: reference for DADA2's `assignTaxonomy()`, headers `>rank1;rank2;...;species`.
    - `*.n<N>.general.fna.gz`: the same sequences with their original headers.
    - `*.n<N>.names.txt`: the names of the exported sequences, one per line, as used in the other files of the run.
  - With `--run_reference`, a reference phylogeny for placing unknown sequences:
    - `*.raxml.bestTree`: the phylogeny, constrained by the taxonomy.
    - `*.raxml.bestModel`: the RAxML-NG model of the tree.
    - `*.reference.fasta`: the alignment the tree was built from.
    - `*.reference.tax`: the declared taxonomy of the same sequences.
    - `*.guide.nwk`: the taxonomy-only guide tree used as the constraint.

</details>

`export/` holds everything the pipeline produces for use after the run.
The exported sequences are the ones that passed every check that ran, at most N per species, ranked by weight times length; with `all` instead of a number, the files are named `*.all.*` and hold every sequence that passed.
See [Reference export](usage.md#reference-export) and [Reference for phylogenetic placement](usage.md#reference-for-phylogenetic-placement) in the usage documentation for how the sequences are chosen and how the reference phylogeny is built.

## Input preparation

### Sequence selection

<details markdown="1">
<summary>Output files</summary>

- `seqkit/`
  - `*.fa`, `*.fa.gz`: the `--sequences` records whose header matches the `--seqgrep` pattern.

</details>

Only with `--seqgrep`. [SeqKit](https://bioinf.shenwei.me/seqkit/) `grep` keeps the records whose whole header line matches the pattern, for example `d__Archaea` to keep only the archaea of a GTDB file. All later steps see only these records.

### Taxonomy resolution

<details markdown="1">
<summary>Output files</summary>

- `resolvetaxonomy/`
  - `*.resolved.tax`: the taxonomy of every sequence.
  - `*.resolved.<ext>`: the sequences in their input format (FASTA, Clustal or PHYLIP), with headers cut down to the sequence name.
  - `*.warnings.txt`: one line per warning, for example when taxonomy text in the headers was ignored because `--taxonomy` was given. Empty if there were none.

</details>

The taxonomy comes from the `--taxonomy` file when one is given.
Without it, it is read from the sequence headers, as in GTDB files (`>id taxonomy;string`).
Empty or placeholder ranks at the end of a lineage (such as `s__`, `;;` or a name matching `--unresolved_rank_pattern`) are removed, so a sequence annotated only to genus ends at the genus.
Sequences without any annotation are dropped, with a warning.

### Name check

<details markdown="1">
<summary>Output files</summary>

- `checknameconsistency/`
  - `*.checked.tax`: the taxonomy, with safe sequence names.
  - `*.checked.<ext>`: the sequences in their input format, with the same safe names.

</details>

This step checks that the taxonomy and the sequences list exactly the same sequence names, and stops the pipeline with a list of the names if they do not.
It then replaces characters that some tools cannot handle (any character other than letters, digits, `_`, `.`, `-`, `|` and `/`) with `_`.
These are the names used in all later files, except the `*.general.fna.gz` files in [`export/`](#export), which keep the original headers.

### Conversion to FASTA

<details markdown="1">
<summary>Output files</summary>

- `emboss/`
  - `*.fasta`: the checked sequences in FASTA format.

</details>

EMBOSS `seqret` converts the sequences to FASTA, the format all later steps read.

## Weighted clustering

The steps in this group reduce the input to one representative per cluster of similar sequences and per taxon before the slower steps run.
They are skipped with `--skip_clustering`.
See [Weighted clustering](usage/weighted_clustering.md) for the design.

### GTDB weights

<details markdown="1">
<summary>Output files</summary>

- `gtdbweighttranslate/`
  - `*.weights.tsv`: a weight for each sequence, tab-separated without a header (sequence name, weight), in the format of `--sequence_weights`.

</details>

Only with `--upstream GTDB`. The weights are computed from the GTDB genome metadata (`--gtdb_bac120_metadata`, `--gtdb_ar53_metadata`): genome category, CheckM quality, type material, contig length and agreement with an independent SILVA classification.
Sequences that do not match a metadata row get no weight here and count as weight 1 later.
An explicit `--sequence_weights` file wins over these weights.

### Weight filter

<details markdown="1">
<summary>Output files</summary>

- `weightfilter/`
  - `*.weightfiltered.tax`: the taxonomy of the sequences that were kept.
  - `*.weightfiltered.fasta`: the sequences that were kept, unchanged.
  - `*.weightfiltered.ungapped.fasta`: the same sequences without gap characters, for VSEARCH.
  - `*.excluded.tsv`: the sequences dropped for a weight below `--min_weight`.

</details>

A sequence is dropped when its weight (`--sequence_weights` or `--upstream GTDB`; 1 for a sequence without a weight) is below `--min_weight`. Without `--min_weight` every sequence is kept, and `*.excluded.tsv` has only its header.
`*.excluded.tsv` is tab-separated, with a header line and these columns:

| Column                 | Description                       |
| ---------------------- | --------------------------------- |
| `seq_name`             | Sequence name                     |
| `weight`               | Weight of the sequence            |
| `min_weight_threshold` | The `--min_weight` that it missed |

### Clustering

<details markdown="1">
<summary>Output files</summary>

- `vsearch/`
  - `*.uc.tsv.gz`: the clusters found by [VSEARCH](https://github.com/torognes/vsearch), in VSEARCH's UC format.

</details>

VSEARCH clusters the sequences at `--cluster_identity` (default 1.0, so only identical sequences end up in the same cluster).
The UC file is tab-separated, without a header, one row per sequence plus one summary row per cluster.
The first column is the row type: `S` for the first sequence of a cluster (the seed), `H` for a sequence added to a cluster, and `C` for the cluster summary.
The second column is the cluster number, the third the sequence length, the fourth the identity to the seed in percent, the ninth the sequence name and the tenth the name of the seed.
See the [VSEARCH manual](https://github.com/torognes/vsearch) for all columns.

### Cluster selection

<details markdown="1">
<summary>Output files</summary>

- `clusterselect/`
  - `*.representatives.fasta`: one sequence per pair of cluster and taxon, the representatives.
  - `*.representatives.tax`: the taxonomy of the representatives.
  - `*.selection.tsv`: every sequence that was considered, its score and whether it was chosen.

</details>

Within each cluster, the sequences are grouped by their declared taxonomy, and the sequence with the highest score in each group becomes the representative.
Only the representatives go on to the next steps. The other members of the group get the same verdict as their representative in the [export](#export).
`*.selection.tsv` is tab-separated, with a header line and these columns:

| Column              | Description                                                                                        |
| ------------------- | -------------------------------------------------------------------------------------------------- |
| `cluster_id`        | VSEARCH cluster number                                                                             |
| `taxon`             | Declared taxonomy, as a `;`-separated lineage                                                      |
| `seq_name`          | Sequence name                                                                                      |
| `length`            | Sequence length, without gaps                                                                      |
| `weight`            | Weight of the sequence, 1 when it has none                                                         |
| `score`             | Length divided by the length of the longest sequence with the same cluster and taxon, times weight |
| `is_representative` | `True` for the sequence chosen as representative, `False` for the others                           |

## raxtax prefilter

The steps in this group classify the sequences with [raxtax](https://pubmed.ncbi.nlm.nih.gov/41259660/), a fast k-mer based classifier, and flag the clearly mislabelled ones before the slower placement step.
They are skipped with `--skip_raxtax`.

### raxtax input

<details markdown="1">
<summary>Output files</summary>

- `raxtaxformat/`
  - `*.raxtax.fasta`: the sequences without gaps, with headers in raxtax's format, `>id;tax=<lineage>;`.

</details>

The same file is used as both the database and the queries of the raxtax run.

### raxtax classification

<details markdown="1">
<summary>Output files</summary>

- `raxtax/user-alignment/`
  - `raxtax.out`: raxtax's classification of every sequence, with a confidence for each rank.
  - `raxtax.log`: the raxtax log.

</details>

Every sequence is classified against all the others.
raxtax runs with `--skip-exact-matches`, so that a sequence does not simply match itself and get its own, possibly wrong, label back.

### raxtax filter

<details markdown="1">
<summary>Output files</summary>

- `raxtaxfilter/`
  - `*.mislabels.tsv`: the sequences flagged by raxtax, with the same columns as the [mislabels table](#mislabels).
  - `*.filtered.fasta`: the sequences that were not flagged.
  - `*.filtered.tax`: the taxonomy of the sequences that were not flagged.

</details>

A sequence is flagged when its best raxtax hit disagrees with its declared taxonomy at the checked rank with a confidence of at least `--raxtax_min_confidence` (default 0.9).
The checked rank is set by `--raxtax_filter_rank`, counted from the end of the lineage (default 2, the genus in a lineage that ends at species).
A sequence is not checked if its lineage is shorter than that rank, or if no other sequence has the same declared taxon at that rank.
Flagged sequences are reported in the [mislabels table](#mislabels) and do not go on to the placement step.

## Alignment

These steps run for every input, but only unaligned sequences are aligned.
See [Aligning the sequences](usage/aligning_sequences.md) for the options.

### Alignment check

<details markdown="1">
<summary>Output files</summary>

- `checkaligned/`
  - `*.aligned.fasta` or `*.unaligned.fasta`: the sequences, unchanged.

</details>

The file name shows what the pipeline found: `aligned` if every sequence has the same length, `unaligned` if not.
Unaligned sequences are aligned in the next step; aligned ones go straight to the [gap filter](#gap-filter).

### Alignment

<details markdown="1">
<summary>Output files</summary>

- `alignment/`
  - `user-alignment.afa`: the aligned sequences, in aligned FASTA format.

</details>

Only for unaligned input.
With `--hmm`, the sequences are aligned to the HMM profile with [HMMER](http://hmmer.org/) `hmmalign`, in chunks of `--align_chunk_size` sequences that are joined afterwards, and the alignment is cut down to the columns of the profile.
Without `--hmm`, all sequences are aligned in one job with [MAFFT](https://mafft.cbrc.jp/alignment/software/).

### HMM profile

<details markdown="1">
<summary>Output files</summary>

- `hmmer/`
  - `hmm.hmm`: the profile named by `--hmm_name`, taken from the `--hmm` file.

</details>

Only with `--hmm_name` and unaligned input. HMMER `hmmfetch` picks the profile out of a file with several profiles.

## Short-sequence filters

Sequences that are too short or incomplete cannot be placed reliably, so one of these filters drops them before the placement step.
Which one runs depends on how the alignment was made.

### Gap filter

<details markdown="1">
<summary>Output files</summary>

- `gapfilter/`
  - `*.gapfiltered.fasta`: the alignment of the sequences that were kept.
  - `*.gapfiltered.tax`: the taxonomy of the sequences that were kept.
  - `*.excluded.tsv`: the sequences that were dropped.

</details>

For an alignment that was given as input or made with MAFFT. Sequences whose share of non-gap positions is below `--min_nongap` (default 0.3) are dropped.
Skipped with `--skip_gapfilter`.
`*.excluded.tsv` is tab-separated, with a header line and these columns:

| Column                 | Description                                       |
| ---------------------- | ------------------------------------------------- |
| `seq_name`             | Sequence name                                     |
| `non_gap_proportion`   | Share of the alignment columns that are not a gap |
| `min_nongap_threshold` | The `--min_nongap` that it missed                 |

### Profile coverage filter

<details markdown="1">
<summary>Output files</summary>

- `profilecover/`
  - `*.covfiltered.fasta`: the alignment of the sequences that were kept, including rescued ones.
  - `*.covfiltered.tax`: the taxonomy of the sequences that were kept.
  - `*.cov_excluded.tsv`: every sequence below `--min_profile_cover`, dropped or rescued.

</details>

For an alignment made with `hmmalign`. Sequences that cover less than `--min_profile_cover` (default 0.8) of the HMM profile are dropped, except that a taxon with no sequence above that value keeps its best-covered sequence, if it reaches `--min_profile_cover_rescue` (default 0.5).
Skipped with `--skip_profile_cover`.
See [Profile coverage and fragment rescue](usage/profile_coverage.md) for the reasons.
`*.cov_excluded.tsv` is tab-separated, with a header line and these columns:

| Column                   | Description                                                                                                |
| ------------------------ | ---------------------------------------------------------------------------------------------------------- |
| `seq_name`               | Sequence name                                                                                              |
| `profile_coverage`       | Share of the profile's columns that the sequence covers                                                    |
| `min_coverage_threshold` | The `--min_profile_cover` that it missed                                                                   |
| `status`                 | `excluded` if the sequence was dropped, `rescued` if it was kept as the best-covered sequence of its taxon |

## Reference tree and placement

The steps in this group build a reference tree and test every sequence against it with the SATIVA method: each sequence is taken out of the tree and placed back, and a sequence that lands outside its declared taxon is flagged.
They are skipped with `--skip_sativa`.

### Taxonomy tree

<details markdown="1">
<summary>Output files</summary>

- `taxonomytree/`
  - `*.guide.nwk`: the taxonomy written as a tree in Newick format, one branch per taxon.

</details>

This tree is the constraint of the RAxML-NG search: the reference tree must keep every taxon together.

### Reference tree

<details markdown="1">
<summary>Output files</summary>

- `raxmlng/`
  - `*.raxml.bestTree`: the reference tree, constrained by the taxonomy tree.
  - `*.raxml.bestModel`: the substitution model and its parameters.
  - `*.raxml.log`: the RAxML-NG log.
  - `*.raxml.mlTrees`, `*.raxml.startTree`: all trees found in the search, and the trees it started from.

</details>

The tree is built by [RAxML-NG](https://github.com/amkozlov/raxml-ng), with the model `--raxmlng_model` (model selection by default) and the options `--raxmlng_fast`, `--raxmlng_seed` and `--raxmlng_blopt`.
With `--run_divided_tree_search`, the tree is searched in parts that are joined into one tree, and the files are from the final RAxML-NG evaluation of that tree; there are no `*.mlTrees` and `*.startTree` files then.
See [Divided tree search](usage/divided_tree_search.md).

### Leave-one-out placement

<details markdown="1">
<summary>Output files</summary>

- `sativaepang/`
  - `*.refjson`: the SATIVA reference: the tree, the alignment and the taxonomy in one file.
  - `*.model`: the model used for the placements.
  - `*.mis`: sativa-epang's own report of mislabelled sequences.
  - `*.mislabels.tsv`: the same mislabels with the columns of the [mislabels table](#mislabels).
  - `*.summary.txt`: run statistics: the number of sequences scored and of putative mislabels.
  - `lootasks/*.l1o_tasks/`: the leave-one-out folds before placement, one `fold_NNN/` directory per fold with the reduced reference (`ref.fasta`, `ref.nwk`), the sequences taken out (`query.fasta`) and the model, plus a `manifest.json` that lists the folds and a `model` file.
  - `looplace/*.l1o_tasks/`: the same folds after placement by [EPA-ng](https://github.com/pierrebarbera/epa-ng), with the placements (`epa_result.jplace`) and logs in each fold directory. Only when all folds are placed in one job.
  - `mergefolds/*.l1o_tasks/`: the placed folds gathered from all jobs. Only when `--folds_per_job` splits the placement over more than one job.

</details>

The placement uses the [sativa-epang](https://github.com/Aaramis/sativa-epang) fork of SATIVA, which places with EPA-ng.
The sequences are split into folds; for each fold, its sequences are removed from the reference tree and placed back, and the placements are scored against the declared taxonomy.

## Reporting

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing the data of the report.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising the run.
For this pipeline the report holds the software versions of every tool, for future traceability, a summary of the parameters and a description of the methods.
For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report_<date>.html`, `execution_timeline_<date>.html`, `execution_trace_<date>.txt` and `pipeline_dag_<date>.html`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `nf_core_taxmarker_software_mqc_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameters are used when running the pipeline.
  - Parameters used by the pipeline run: `params_<date>.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
