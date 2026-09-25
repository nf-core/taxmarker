# nf-core/taxmarker: Weighted clustering

Before anything else sees the input (alignment, raxtax, placement), the pipeline reduces it to one representative sequence per (cluster, taxon) pair.
This collapses near-duplicate sequences of the same taxon -- common in real reference databases such as GTDB, which can carry dozens of genomes per species -- before the much more expensive downstream steps run on them.

Three parameters control this, all optional:

- `--sequence_weights`: a two-column table (`sequence_id<TAB>weight`) used to break ties when picking a representative.
  Sequences absent from the table default to weight 1.
  Omit entirely and every sequence is weight 1, which degrades to "pick the longest sequence per (cluster, taxon) pair".
  The clustering logic itself is generic and doesn't care where the weights came from -- see `--upstream` below for a way to compute it automatically instead of supplying it directly.
- `--min_weight`: an absolute cutoff applied before clustering.
  Sequences below it are excluded outright, regardless of `--sequence_weights`.
  Has no effect unless both this and `--sequence_weights` are set -- relative per-cluster selection alone can't exclude a badly-mislabeled sequence that never co-clusters with anything (e.g. a singleton cluster).
- `--cluster_identity`: the [VSEARCH](https://github.com/torognes/vsearch) clustering identity threshold, 0-1, default `1.0`.
  The default is a safe no-op: pure dereplication of exact/near-identical duplicates.
  Lower it for more aggressive input-size reduction on very large datasets.

## How VSEARCH decides what counts as "identical"

`--cluster_identity` doesn't mean what a quick reading suggests, so it's worth spelling out exactly what VSEARCH does, confirmed against its own manual (`man vsearch`, v2.32.0):

- VSEARCH always aligns the full length of both sequences against each other with a global Needleman-Wunsch algorithm (full dynamic programming) -- never a local/Smith-Waterman alignment restricted to the best-matching region.
- Terminal gaps (a dangling, non-overlapping end on either sequence) get a much cheaper gap penalty than internal gaps by default (open 2 vs. 20, extend 1 vs. 2), so a non-overlapping end is simply gapped out rather than forced into mismatches.
- The identity score itself (`--iddef 2`, VSEARCH's own default) is `matches / (alignment length - terminal gaps)`: terminal-gap columns are excluded from **both** the numerator and the denominator.

The consequence: identity reflects only the _overlapping_ region between two sequences, and the overlap's length relative to either full sequence plays no role in the score.
Two sequences that only partially overlap -- each with its own unique, non-overlapping end -- can still reach 100% identity and be clustered together, provided the shared middle matches perfectly, no matter how small a fraction of either sequence's total length that shared region is.
Conversely, mismatches that fall _inside_ the shared/aligned region (not at the ends) always count against identity, so genuinely divergent sequences are correctly kept apart at the pipeline's default `--cluster_identity 1.0`.

This is what lets the pipeline's default settings correctly merge, say, a full-length 16S sequence and a shorter deposited fragment of the same organism (the fragment's implied "missing" tail is just a terminal gap on the longer sequence).
The corresponding risk -- a short, spuriously-identical overlap merging two sequences that shouldn't be merged -- is bounded downstream: representative selection groups by `(cluster, taxon)`, so a cross-taxon spurious merge can't produce a wrong cross-taxon representative; only a same-taxon spurious merge could pick between two genuinely different same-species fragments, which the pipeline treats as an acceptable simplification for its purpose (one representative sequence per taxon, not preservation of every individual record).

## `--upstream`: computing weights automatically for a known data source

`--upstream` names the data source `--sequences`/`--taxonomy` came from, to enable source-specific pipeline behaviour.
Currently the only recognised value is `GTDB`.

### `--upstream GTDB`

Requires `--gtdb_bac120_metadata` and/or `--gtdb_ar53_metadata` -- GTDB's own per-release `bac120_metadata_*.tsv[.gz]`/`ar53_metadata_*.tsv[.gz]` files, unmodified.
Only `--sequences` records matching GTDB's own `ssu_all` naming convention (`<accession>~<contig>`, e.g. `RS_GCF_002194975.1~NZ_NHWV01000143.1`) with a matching metadata row get a weight; everything else falls back to the generic default of weight 1.

If `--sequence_weights` is also given, that file wins (with a warning) -- the same precedent as an explicit `--taxonomy` file winning over taxonomy embedded in `--sequences` headers.

The formula, converged on [issue #15](https://github.com/nf-core/taxmarker/issues/15) and calibrated against real GTDB r226 metadata:

```
weight = category_term + quality_score + type_material_term + contig_len_term + mismatch_penalty

quality_score:    checkm2_completeness - 5*checkm2_contamination (falls back to checkm_* if checkm2_* is missing)
category_term:    isolate (ncbi_genome_category = none) = 10, SAG (derived from single cell) = 5,
                   MAG (derived from metagenome) = 0, other/unknown = 0
type_material_term (ncbi_type_material_designation): na = 0, type material = 5,
                   synonym type material = 3, pathotype = 3, neotype/reftype/clade exemplar = 2
contig_len_term:   this sequence's own header `[contig_len=...]` -- <2000bp = 0, 2000-20000bp = 1,
                   20000-200000bp = 2, >=200000bp = 3, missing = 0
mismatch_penalty:  shallowest rank (domain->genus) where gtdb_taxonomy disagrees with a
                   synonym-normalised ssu_silva_taxonomy: domain = -200, phylum = -150,
                   class = -40, order = -25, family = -10, genus = -5, no evidence = 0
```

The `mismatch_penalty` term catches a genome whose declared (GTDB, marker-protein-based) taxonomy disagrees with an independent classification of its own extracted 16S sequence (SILVA, BLAST-based) -- a strong, cheap signal for a contaminating or mislabeled 16S copy that checkm2's own quality score is blind to (it's marker-protein based, so a perfect-looking genome can still carry a bad 16S).
The penalty magnitude at each rank is deliberately larger than the maximum plausible positive total from the other four terms (118), so a genome with a shallow-rank mismatch reliably drops below `--min_weight 0` regardless of how good it otherwise looks -- **using `--upstream GTDB` without also setting `--min_weight 0` means these mismatches are computed but never actually excluded.**

GTDB and SILVA independently revise taxonomic nomenclature over time (e.g. `Desulfobacterota`/`Thermodesulfobacteriota` are the same phylum under different names), so a naive name comparison would double-count nomenclature drift as contamination signal.
To avoid this, synonym sets are built per rank from isolate-only GTDB/SILVA co-occurrence in the same metadata file: a SILVA name occurring for at least 5% of a GTDB name's isolate rows at a given rank is treated as an accepted synonym there, checked empirically against real GTDB r226 data (true synonyms sit at ~90-100% co-occurrence, genuine errors under ~1%).
