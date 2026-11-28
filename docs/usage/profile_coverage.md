# nf-core/taxmarker: Profile coverage and fragment rescue

When `--sequences` is unaligned, the pipeline aligns it to an HMM profile (`--hmm`) and masks the alignment down to the profile's match-state columns.
The proportion of those columns a sequence covers is its _profile coverage_: 1.0 for a complete gene, lower for a partial one.

Two parameters control what happens to partial sequences:

- `--min_profile_cover` (default `0.8`): sequences at or above this coverage are kept.
- `--min_profile_cover_rescue` (default `0.5`): if a taxon has no sequence at or above `--min_profile_cover`, its best-covered sequence is kept anyway, provided it reaches this floor.
  Only one sequence per taxon is rescued, and a taxon is the full taxonomy string, i.e. usually a species.
  Set it to `--min_profile_cover` or higher to turn the rescue off.

`profilecover/*.cov_excluded.tsv` lists every sequence below `--min_profile_cover`, with its coverage and a `status` of `excluded` or `rescued`.
`--skip_profile_cover` turns the filter off altogether.

## Why rescue fragments

Many taxa, especially those known only from metagenome-assembled genomes, have only partial marker-gene sequences.
Dropping them all removes those taxa from the reference entirely, and a fragment still carries phylogenetic signal.

Measured on the archaeal 16S sequences of GTDB r226, after weighted clustering at the default settings (8863 representatives):

- Coverage of the partial sequences is spread evenly from 0.1 to 0.8, with no natural cut-off between complete and partial sequences.
- `--min_profile_cover 0.8` alone keeps 4583 sequences, and loses every sequence of 36% of species (3044 of 4774 kept) and 27% of genera (1242 of 1698).
  Lowering the threshold alone does not solve this: at 0.5, 23% of species are still lost.
- Rescue at the default floor of 0.5 adds 643 sequences, and keeps 3687 species, 1461 genera and 473 of 523 families.

## Do fragments carry consistent signal?

Three checks on the same data, rescuing down to a floor of 0.3 (1100 fragments) to see how coverage matters:

- **Unconstrained tree.**
  In a FastTree (GTR+Γ) tree built with and without the fragments, the share of a tip's sister clade that belongs to the same genus averages 0.90 for complete sequences.
  It is 0.81 for fragments with coverage 0.5-0.8 and 0.73 for those at 0.3-0.5, which is why the default floor is 0.5.
  At family level the figures are 0.94, 0.89 and 0.87.
  Adding the fragments costs 37 of 597 genera their monophyly (7 gain it), and gives 156 more genera at least two sequences.
- **Placement into a tree of complete sequences.**
  Placed with EPA-ng into the pipeline's own reference tree built without them, 87% of fragments with coverage 0.5-0.8 land inside or at the stem of their declared genus (79% at 0.3-0.5), and 95% inside their family.
  The median like-weight ratio of the best placement is 1.0.
- **Branch lengths.**
  Fragments have a median terminal branch of 0.017, against 0.003 for complete sequences.
  Missing data does not cause this: giving 1100 complete sequences the gap patterns of the fragments leaves their terminal branches unchanged.
  Most of the difference comes from sampling instead: a rescued fragment is by construction the only sequence of its species, and complete sequences that are alone in their species have a median of 0.011.

A small share of fragments disagree with their declared genus, more than complete sequences do.
Some of that is real mislabelling rather than noise, and the leave-one-out placement is there to catch it.

## Interaction with the raxtax prefilter

The raxtax prefilter runs before alignment, so it classifies partial sequences before this filter sees them.
Its flags on sequences that would later have been dropped for low coverage still count as mislabels.
