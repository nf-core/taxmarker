/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENSURE_ALIGNED - transparently align unaligned input via hmmalign or MAFFT

    Sequences that are already a multiple sequence alignment (every record the same
    length) pass through unchanged. Sequences that aren't are aligned against an HMM
    profile (hmm, optionally a specific named profile within it via hmm_name) using
    hmmalign, then converted back to plain FASTA. Without an HMM they are aligned
    from scratch with MAFFT instead. No separate mode-switch param needed,
    detection is by content.

    Workflow:
      1. Detect aligned vs unaligned (equal sequence lengths?)     (CHECKALIGNED)
      2. unaligned, no hmm: align all sequences in one job         (MAFFT_ALIGN)
      Steps 3-8 are for unaligned input with hmm:
      3. extract the named profile, if given                       (HMMER_HMMFETCH)
      4. split into chunks of chunk_size sequences                 (splitFasta)
      5. per chunk: align against the profile                      (HMMER_HMMALIGN)
      6. per chunk: mask down to the profile's own match-state
         columns, discarding insert-state columns                  (HMMER_ESLALIMASK)
      7. per chunk: convert Stockholm back to aligned FASTA        (HMMER_ESLREFORMAT)
      8. concatenate the chunks in input order and
         decompress (eslreformat always gzips)                     (FIND_CONCATENATE)

    Chunking bounds hmmalign's memory, which grows with the number of sequences
    (about 60 kB each) because it holds every trace until the alignment is written.
    Every chunk is masked to the same match-state columns, so the masked chunks
    concatenate into exactly the alignment a single hmmalign run would give.

    Masking matters beyond tidiness: hmmalign's raw output width is match-state
    columns *plus* the union of every sequence's own insert-state columns, so a
    single divergent sequence's insertion pads every other sequence with extra
    gap columns there. Left unmasked, that badly deflates a downstream non-gap-
    proportion filter for genuinely full-length sequences too (confirmed
    empirically: masking took the median non-gap proportion across a real
    121-sequence archaeal 16S set from 0.61 to 0.99). Masking first makes the
    alignment width exactly the HMM's match-state count.

    MAFFT aligns all sequences in one job, so its time and memory grow with the
    input: it suits up to some thousands of sequences. Larger inputs need an HMM.

    The alignments are emitted SEPARATELY (alignment_native, alignment_from_hmm),
    not merged into one, because the caller applies a different gap/coverage filter
    -- with its own threshold default -- to each: an alignment whose columns come
    from the sequences themselves (given, or made by MAFFT) is naturally wide with
    real cross-lineage indels even for full-length sequences, whereas this
    subworkflow's masked hmmalign output has a fixed-length, much tighter
    distribution (see GAPFILTER vs PROFILECOVER in workflows/taxmarker.nf).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECKALIGNED      } from '../../../modules/local/checkaligned/main'
include { MAFFT_ALIGN       } from '../../../modules/nf-core/mafft/align/main'
include { HMMER_HMMFETCH    } from '../../../modules/nf-core/hmmer/hmmfetch/main'
include { HMMER_HMMALIGN    } from '../../../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_ESLALIMASK  } from '../../../modules/nf-core/hmmer/eslalimask/main'
include { HMMER_ESLREFORMAT } from '../../../modules/nf-core/hmmer/eslreformat/main'
include { FIND_CONCATENATE  } from '../../../modules/nf-core/find/concatenate/main'

workflow ENSURE_ALIGNED {

    take:
    ch_sequences // channel: sequences file, aligned or not, already FASTA (format-normalised by the caller)
    hmm          // value:   path to an HMM profile database, or null/empty to align unaligned input with MAFFT
    hmm_name     // value:   name of a specific profile within hmm, or null/empty
    chunk_size   // value:   sequences aligned per hmmalign job

    main:
    CHECKALIGNED(ch_sequences.map { [ [ id: 'user-alignment' ], it ] })

    // hmm is a plain value, known before execution starts, so unaligned input can be
    // routed at compose time: the other branch gets an empty channel and never runs.
    def ch_unaligned_to_hmm   = hmm ? CHECKALIGNED.out.unaligned : channel.empty()
    def ch_unaligned_to_mafft = hmm ? channel.empty() : CHECKALIGNED.out.unaligned
        .map { meta, fasta ->
            def n = fasta.countFasta()
            if (n > 20000) {
                log.warn("Aligning ${n} sequences with MAFFT in one job, which may take very long and a lot of memory. Give --hmm to align in chunks.")
            }
            [ meta, fasta ]
        }

    if (hmm_name && !hmm) {
        log.warn("--hmm_name is ignored without --hmm.")
    }

    // compress=false: the module's gzip pipe would hide a failed MAFFT behind pigz's exit status.
    MAFFT_ALIGN(
        ch_unaligned_to_mafft,
        [ [], [] ], [ [], [] ], [ [], [] ], [ [], [] ], [ [], [] ],
        false
    )

    def ch_hmm_profile
    if (hmm) {
        def ch_hmm_db = channel.fromPath(hmm, checkIfExists: true)
            .map { [ [ id: 'hmm' ], it ] }
        if (hmm_name) {
            // hmmfetch works directly against a raw (unindexed) multi-profile
            // database when given an explicit key -- no separate --index step needed.
            HMMER_HMMFETCH(ch_hmm_db, hmm_name, [], [])
            ch_hmm_profile = HMMER_HMMFETCH.out.hmm
        } else {
            ch_hmm_profile = ch_hmm_db
        }
    } else {
        ch_hmm_profile = channel.empty()
    }

    // Zero-padded chunk number in the id: FIND_CONCATENATE joins its inputs in file name order.
    def ch_chunks = ch_unaligned_to_hmm
        .map { _meta, fasta -> fasta }
        .splitFasta(by: chunk_size, file: true)
        .toList()
        .flatMap { chunks ->
            chunks.withIndex().collect { chunk, i -> [ [ id: "user-alignment.${String.format('%06d', i)}" ], chunk ] }
        }

    // first(): a value channel, so every chunk pairs with the one profile.
    HMMER_HMMALIGN(ch_chunks, ch_hmm_profile.map { _meta, hmm_file -> hmm_file }.first())

    // --rf-is-mask derives the mask from the alignment's own #=GC RF (reference/
    // match-state) annotation -- no separate maskfile needed, unlike the tool's
    // other usage form. The six false values are the module's optional
    // mask-report-file toggles (fmask/gmask/pmask), none of which are needed here.
    HMMER_ESLALIMASK(
        HMMER_HMMALIGN.out.sto.map { meta, sto -> [ meta, sto, false, false, false, false, false, false ] },
        []
    )

    HMMER_ESLREFORMAT(HMMER_ESLALIMASK.out.maskedaln, '')

    // eslreformat always gzips; the .gz-less prefix (conf/modules.config) makes
    // FIND_CONCATENATE decompress, matching the plain-text native branch.
    FIND_CONCATENATE(
        HMMER_ESLREFORMAT.out.seqreformated
            .toList()
            .filter { chunks -> chunks }
            .map { chunks -> [ [ id: 'user-alignment' ], chunks.collect { _meta, fasta -> fasta } ] }
    )

    emit:
    alignment_native      = CHECKALIGNED.out.aligned.mix(MAFFT_ALIGN.out.fas).map { _meta, fasta -> fasta } // channel: alignment file, already aligned on input or aligned by MAFFT
    alignment_from_hmm    = FIND_CONCATENATE.out.file_out.map { _meta, fasta -> fasta } // channel: alignment file, aligned via hmmalign then masked to match-state columns
}
