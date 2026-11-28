/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RAXTAX_PREFILTER - optional fast prefilter ahead of alignment and the expensive
    EPA-ng-based SATIVA subworkflow, using raxtax (https://github.com/noahares/raxtax)
    to self-classify the reference set and catch severely mislabeled sequences early.

    Runs on unaligned sequences: raxtax classifies plain sequences and never needs an
    alignment, so running this before ENSURE_ALIGNED lets flagged sequences skip
    alignment entirely, not just placement.

    Workflow:
      1. Rewrite headers to raxtax's format (also strips any
         gap characters, a no-op on unaligned input)         (RAXTAXFORMAT)
      2. Self-classify (database == query)                   (RAXTAX)
      3. Keep the best hit per query, flag disagreements
         at params.raxtax_filter_rank                         (RAXTAXFILTER)

    Sequences are expected to already be normalised to FASTA by the caller
    (workflows/taxmarker.nf) -- this subworkflow no longer does that itself.

    Sequences RAXTAXFILTER flags never reach alignment or EPA-ng placement: they're
    reported directly in the final mislabels output (see workflows/taxmarker.nf), tagged
    method=raxtax to distinguish them from SATIVASCORE's own (method=sativa) rows.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RAXTAXFORMAT  } from '../../../modules/local/raxtaxformat/main'
include { RAXTAX        } from '../../../modules/local/raxtax/main'
include { RAXTAXFILTER  } from '../../../modules/local/raxtaxfilter/main'

workflow RAXTAX_PREFILTER {

    take:
    ch_taxonomy  // channel: taxonomy file (seq_name<TAB>rank1;rank2;...)
    ch_sequences // channel: unaligned sequences file, already normalised to FASTA by the caller

    main:
    // Give both inputs a shared meta so they can be joined back together below;
    // matches the fixed 'user-alignment' id used throughout subworkflows/local/sativa.
    def ch_meta_taxonomy    = ch_taxonomy.map  { [ [ id: 'user-alignment' ], it ] }
    def ch_meta_sequences   = ch_sequences.map { [ [ id: 'user-alignment' ], it ] }

    // Rewrite headers to raxtax's `>id;tax=<lineage>;` form.
    RAXTAXFORMAT(ch_meta_sequences.join(ch_meta_taxonomy))
    def ch_raxtax_fasta = RAXTAXFORMAT.out.fasta

    // Self-classification: database and query are the same file. --skip-exact-matches
    // (set via conf/modules.config) drops each query's own trivial self-hit, so the
    // classification reflects genuine similarity to the *other* reference sequences.
    RAXTAX(ch_raxtax_fasta, ch_raxtax_fasta.map { _meta, fasta -> fasta })

    RAXTAXFILTER(
        RAXTAX.out.out
            .join(ch_meta_taxonomy)
            .join(ch_meta_sequences)
    )

    emit:
    taxonomy  = RAXTAXFILTER.out.taxonomy.map  { _meta, tax -> tax }  // channel: taxonomy file, sequences NOT flagged by raxtax
    sequences = RAXTAXFILTER.out.alignment.map { _meta, seq -> seq }  // channel: sequences file, sequences NOT flagged by raxtax
    mislabels = RAXTAXFILTER.out.mislabels                            // [ meta, tsv ]  raxtax-flagged sequences, method=raxtax
}
