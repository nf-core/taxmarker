/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WEIGHTED_CLUSTERING - reduce a large input sequence set before alignment/placement

    Design: nf-core/taxmarker#15. A generic per-sequence weight (sequence_id<TAB>weight,
    defaulting to 1 for every sequence if --sequence_weights isn't given -- the
    multiplicative identity, degrading gracefully to "pick the longest sequence per
    taxon") plus VSEARCH clustering by similarity.

    Workflow:
      1. Strip gap characters (caller doesn't guarantee unaligned input --
         ENSURE_ALIGNED, which detects aligned-vs-not, runs later), then
         apply the min_weight cutoff -- absolute exclusion, before
         clustering (saves clustering compute on sequences that would be
         dropped anyway). Cutoff is a no-op if --min_weight isn't set. (WEIGHTFILTER)
      2. Cluster by similarity, at --cluster_identity (default 1.0,
         pure dereplication -- lower it for more aggressive reduction). (VSEARCH_CLUSTER)
      3. Per (cluster, taxon) pair, pick the representative that
         maximises (length / max-length-in-cluster) x weight.          (CLUSTERSELECT)

    Downstream (raxtax, alignment, placement) sees representatives only -- that's the
    actual point of clustering early. Runs before RAXTAX_PREFILTER, so the raxtax
    prefilter gets the same input-size reduction for free.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WEIGHTFILTER    } from '../../../modules/local/weightfilter/main'
include { VSEARCH_CLUSTER } from '../../../modules/nf-core/vsearch/cluster/main'
include { CLUSTERSELECT   } from '../../../modules/local/clusterselect/main'

workflow WEIGHTED_CLUSTERING {

    take:
    ch_taxonomy         // channel: taxonomy file (seq_name<TAB>rank1;rank2;...)
    ch_sequences        // channel: unaligned sequences file, already normalised to FASTA by the caller
    ch_sequence_weights // channel: single item, itself a 1-element list wrapping the weight-table
                         //          path (or [] if not supplied) -- see the .combine() below for why
    min_weight          // value:   absolute weight cutoff, or null/empty to skip it

    main:
    def ch_meta_taxonomy  = ch_taxonomy.map  { [ [ id: 'user-alignment' ], it ] }
    def ch_meta_sequences = ch_sequences.map { [ [ id: 'user-alignment' ], it ] }

    // WEIGHTFILTER always runs -- it's also where gap characters get stripped, which
    // every sequence needs before VSEARCH_CLUSTER regardless of whether a cutoff is
    // active. -Infinity keeps a real weight >= it always, making the cutoff itself a
    // true no-op when --min_weight isn't set.
    def effective_min_weight = (min_weight != null && min_weight.toString() != '') ? min_weight : Double.NEGATIVE_INFINITY
    // ch_sequence_weights' items are wrapped in an outer 1-element list: .combine()
    // would otherwise splice a bare [] payload's own elements in (zero of them),
    // desyncing this tuple's cardinality instead of giving it one empty-list slot
    // (same pitfall RESOLVETAXONOMY's own optional-taxonomy handling documents) --
    // confirmed empirically that combine() correctly unwraps the outer list back to
    // the plain value (path or []) as it splices it in.
    WEIGHTFILTER(
        ch_meta_taxonomy.join(ch_meta_sequences)
            .combine(ch_sequence_weights)
            .map { meta, tax, seq, weights -> [ meta, tax, seq, weights, effective_min_weight ] }
    )

    // Cluster on the ungapped copy (VSEARCH needs plain sequences); CLUSTERSELECT
    // reads the original (still gapped, if it was) sequences, so already-aligned
    // input still looks aligned to ENSURE_ALIGNED once representatives are picked.
    VSEARCH_CLUSTER(WEIGHTFILTER.out.ungapped_sequences)

    CLUSTERSELECT(
        WEIGHTFILTER.out.taxonomy
            .join(WEIGHTFILTER.out.sequences)
            .join(VSEARCH_CLUSTER.out.uc)
            .combine(ch_sequence_weights)
            .map { meta, tax, seq, uc, weights -> [ meta, tax, seq, uc, weights ] }
    )

    emit:
    taxonomy  = CLUSTERSELECT.out.taxonomy.map  { _meta, tax -> tax } // channel: taxonomy file, representatives only
    sequences = CLUSTERSELECT.out.sequences.map { _meta, seq -> seq } // channel: sequences file, representatives only
}
