# nf-core/taxmarker: Aligning the sequences

The placement test needs a multiple sequence alignment.
The pipeline checks `--sequences` for you: if every record has the same length, it is taken as an alignment, and otherwise it is aligned.
There is no parameter to switch between the cases, only `--hmm`:

| Input                    | What happens                                           | Filter on what is left                   |
| ------------------------ | ------------------------------------------------------ | ---------------------------------------- |
| Aligned                  | Used as it is.                                         | Non-gap proportion (`--min_nongap`)      |
| Unaligned, `--hmm` given | Aligned to the HMM profile with `hmmalign`, in chunks. | Profile coverage (`--min_profile_cover`) |
| Unaligned, no `--hmm`    | Aligned from scratch with MAFFT, in one job.           | Non-gap proportion (`--min_nongap`)      |

In both unaligned cases the alignment is published as `alignment/user-alignment.afa`.
Sequences removed by the filter are listed in `gapfilter/` or `profilecover/`, with the value that excluded them.

## With an HMM profile

Give `--hmm` (and `--hmm_name` if the file holds more than one profile) when a profile exists for your marker.
Barrnap's rRNA profiles are an example.
Each sequence is aligned to the profile on its own, so the result does not depend on which other sequences are in the input.
The alignment is cut down to the profile's match-state columns, which fixes its width and lets the profile coverage of each sequence be measured directly (see [Profile coverage and fragment rescue](profile_coverage.md)).

`hmmalign` needs about 60 kB of memory per sequence, so the sequences are aligned in chunks of `--align_chunk_size` (default 10000) in separate jobs.
The chunks are masked to the same columns and concatenated in input order, so the alignment is the same whatever the chunk size.
This is the choice for large inputs: a run of several hundred thousand 16S sequences was aligned this way.

## Without an HMM profile

When `--hmm` is not given, all unaligned sequences are aligned together with [MAFFT](https://mafft.cbrc.jp/alignment/software/) (`mafft --auto --preservecase`).
`--auto` lets MAFFT choose the method from the number and length of the sequences, and `--preservecase` keeps the case of the input.
This needs no knowledge of the marker, so it also works where there is no profile, for example ITS.

MAFFT aligns everything in one job (`process_high`), so time and memory grow with the number and length of the sequences.
One measurement: 19,990 ITS1 sequences (median length 223 nt) took 10 minutes and 8.5 GB with MAFFT's FFT-NS-2 method.
We have not run MAFFT on inputs of several hundred thousand sequences and do not expect it to be practical there; use an HMM for those.

Unlike `hmmalign`, MAFFT lets the sequences decide the columns, so the alignment can be much wider than a single sequence and most columns hold few residues.
That is why the filter is the same non-gap filter as for a given alignment, with its low default (`--min_nongap 0.3`), and not the profile coverage filter.
The width depends on all the sequences: a few long insertions at different places add many columns that most sequences have as gaps, and can push most sequences below `--min_nongap`.
If the filter removes far more sequences than you expect (`gapfilter/` lists them), lower `--min_nongap` or remove the outliers.
Raise it if your sequences are close relatives and a stricter filter makes sense.

To change MAFFT's options, set `ext.args` of `MAFFT_ALIGN` in a custom configuration file (`-c`), which replaces the default `--auto --preservecase` (without `--preservecase` MAFFT writes the alignment in lower case):

```nextflow
process {
    withName: 'MAFFT_ALIGN' {
        ext.args = '--retree 2 --maxiterate 2 --preservecase'
    }
}
```

### Things to know

- MAFFT does not reverse-complement sequences unless you ask it to, and its `--adjustdirection` renames the sequences it turns around, which would no longer match `--taxonomy`.
  A sequence on the other strand is aligned without error, mostly as gaps and mismatches, and can still pass the gap filter, so make sure all sequences are on the same strand before you run the pipeline.
- Whether the input is aligned is decided by length only: records that are all the same length, or a single record, are taken as an alignment even if they are not.
  Pad or align such input yourself if that is wrong.
- Without `--hmm`, the pipeline warns if the input holds more than 20,000 sequences, since all of them go to one MAFFT job.
- Sequences that are too different to align reliably, such as a marker with very long insertions in some taxa or sequences from distantly related groups, end up in columns shared with few others.
  The test then only has close relatives to compare them with.
  How to handle markers like this, where one global alignment is not meaningful, is discussed in [issue #68](https://github.com/nf-core/taxmarker/issues/68).
- Your own alignment is always an option.
  Aligned input skips this step entirely, so you can use any aligner and settings you prefer and give the result as `--sequences`.
