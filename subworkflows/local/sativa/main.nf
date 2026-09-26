/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SATIVA - Semi-Automatic Taxonomy Improvement and Validation Algorithm
    Algorithm originally from https://github.com/amkozlov/sativa. Builds a
    taxonomy-constrained ML reference tree with RAxML-NG, then delegates
    leave-one-out placement/scoring to Auguste Gardette's sativa-epang fork
    (https://github.com/Aaramis/sativa-epang), which places every held-out
    sequence via EPA-ng.

    Workflow:
      1. Build a taxonomy-constrained ML reference tree     (taxonomy2phylogeny)
      2. Build a sativa-epang reference from that tree      (sativaepang/reference)
      3. Deal the reference into leave-one-out folds        (sativaepang/lootasks)
      4. Place every fold via EPA-ng                        (sativaepang/looplace)
      5. Score placements, report mismatches                (sativaepang/looscore)
      6. Translate the .mis report into this pipeline's own
         mislabels.tsv/summary.txt schema                   (sativaepang/misreport)

main.nf
  └── PIPELINE_INITIALISATION   (subworkflows/local/utils_nfcore_taxmarker_pipeline/main.nf)
  └── NFCORE_TAXMARKER
        └── TAXMARKER            (workflows/taxmarker.nf)  ← main logic lives here
              ├── RESOLVETAXONOMY       (modules/local/resolvetaxonomy/)
              ├── CHECKNAMECONSISTENCY  (modules/local/checknameconsistency/)
              ├── EMBOSS_SEQRET         (modules/nf-core/emboss/seqret/) -- normalises to FASTA
              ├── WEIGHTED_CLUSTERING   (subworkflows/local/weighted_clustering/) -- reduce
              │     the input to one representative per (cluster, taxon) pair before
              │     raxtax/alignment/placement all see it (see nf-core/taxmarker#15)
              ├── RAXTAX_PREFILTER      (subworkflows/local/raxtax_prefilter/) -- optional,
              │     params.skip_raxtax to disable; fast self-classification triage on
              │     unaligned sequences that drops severely mislabeled sequences before
              │     alignment, reporting them directly instead
              ├── ENSURE_ALIGNED        (subworkflows/local/ensure_aligned/) -- transparently
              │     aligns unaligned input via hmmalign (params.hmm); already-aligned
              │     input passes through unchanged
              ├── GAPFILTER / PROFILECOVER (modules/local/{gapfilter,profilecover}/)
              ├── SATIVA (this subworkflow)
              └── MULTIQC         (modules/nf-core/multiqc/)
  └── PIPELINE_COMPLETION        (subworkflows/local/utils_nfcore_taxmarker_pipeline/main.nf)
        sends email / completion summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TAXONOMY2PHYLOGENY      } from '../taxonomy2phylogeny/main'
include { SATIVAEPANG_REFERENCE   } from '../../../modules/local/sativaepang/reference/main'
include { SATIVAEPANG_LOOTASKS    } from '../../../modules/local/sativaepang/lootasks/main'
include { SATIVAEPANG_LOOPLACE    } from '../../../modules/local/sativaepang/looplace/main'
include { SATIVAEPANG_MERGEFOLDS  } from '../../../modules/local/sativaepang/mergefolds/main'
include { SATIVAEPANG_LOOSCORE    } from '../../../modules/local/sativaepang/looscore/main'
include { SATIVAEPANG_MISREPORT   } from '../../../modules/local/sativaepang/misreport/main'

// ─── Subworkflow ──────────────────────────────────────────────────────────────

workflow SATIVA {

    take:
    ch_taxonomy   // channel: [ val(meta), path(taxonomy.tsv) ]
                  //   Tab-separated: seq_name <TAB> Kingdom;Phylum;Class;...
                  //   The taxonomic code (BAC/BOT/ZOO/VIR) is the first token.
                  //   CHECKNAMECONSISTENCY (workflows/sativa.nf) has already verified
                  //   this names the same sequences as ch_alignment and rewritten any
                  //   problematic characters (e.g. parens) by the time it gets here.

    ch_alignment  // channel: [ val(meta), path(alignment) ]
                  //   Aligned, labeled sequences in FASTA format -- normalised to
                  //   FASTA and, if it arrived unaligned, aligned via hmmalign, both
                  //   by the caller (workflows/sativa.nf; see EMBOSS_SEQRET and
                  //   ENSURE_ALIGNED there) before this subworkflow ever sees it.
                  //   Sequence IDs must match the first column of ch_taxonomy.

    taxcode       // value:   sativa-epang taxonomic code: bac/bot/zoo/vir

    folds_per_job // value:   folds each leave-one-out placement job places, or null for all in one

    ch_ref_tree   // channel: [ val(meta), path(tree.nwk) ]
                  //   Pre-built reference tree. Pass Channel.empty() to build one.

    ch_ref_model  // channel: [ val(meta), path(model.txt) ]
                  //   RAxML-NG model file matching ch_ref_tree. Channel.empty() if none.

    main:
    // ch_alignment is already FASTA (normalised once by the caller); give it a meta
    // for the joins/tuples below.
    def ch_alignment_meta = ch_alignment.map { [ [ id: 'user-alignment' ], it ] }

    // ── Phase 1: Reference tree construction ───────────────────────────────────
    //
    // Build a multifurcating guide tree from taxonomy strings, then run RAxML-NG
    // with that tree as a topology constraint and its own automatic model testing
    // (MOOSE, triggered by "DNA" for these nucleotide marker genes). The resulting
    // tree + model are reusable across runs (pass via ch_ref_tree / ch_ref_model to
    // skip this phase).

    // TODO: give externally supplied ch_ref_tree / ch_ref_model precedence over what
    // we just built, once main.nf actually exposes a way to pass them in (currently
    // always called with `[]`, so mixing them in here would inject a spurious
    // empty-list item into the channel).
    def ch_taxonomy_meta = ch_taxonomy.map { [ [ id: 'user-alignment' ], it ] }
    TAXONOMY2PHYLOGENY(
        ch_taxonomy_meta
            .join(ch_alignment_meta)
            .map { meta, taxonomy, alignment -> [ meta, taxonomy, alignment, 'DNA' ] }
    )

    def ch_tree  = TAXONOMY2PHYLOGENY.out.tree
    def ch_model = TAXONOMY2PHYLOGENY.out.model

    // ── Phase 2: sativa-epang reference + leave-one-out scatter ────────────────
    //
    // Hand the RAxML-NG tree+model to sativa-epang via -reftree/-refmodel, skipping
    // its own constrained RAxML search entirely. Joined by meta.id (not paired
    // positionally) before splitting back into sativaepang/reference's three
    // positional inputs, since -reftree/-refmodel carry no meta of their own -- see
    // this project's own Nextflow channel-joining doctrine.
    def ch_reference_input = ch_alignment_meta
        .join(ch_taxonomy_meta)
        .map { meta, alignment, taxonomy -> [ meta, alignment, taxonomy, taxcode ] }
        .join(ch_tree)
        .join(ch_model)
    // ch_reference_input: [ meta, alignment, taxonomy, taxcode, reftree, refmodel ]

    SATIVAEPANG_REFERENCE(
        ch_reference_input.map { meta, alignment, taxonomy, taxcode_item, _reftree, _refmodel -> [ meta, alignment, taxonomy, taxcode_item ] },
        ch_reference_input.map { meta, _alignment, _taxonomy, _taxcode, reftree, _refmodel -> reftree },
        ch_reference_input.map { meta, _alignment, _taxonomy, _taxcode, _reftree, refmodel -> refmodel }
    )

    SATIVAEPANG_LOOTASKS(SATIVAEPANG_REFERENCE.out.refjson)

    // One job per folds_per_job folds, sliced from the fold count lootasks actually
    // wrote into its manifest: deriving the job count from the work means no job is
    // ever handed an empty range. shard and nshards reach the tool as -folds via
    // conf/modules.config, since modules may not read custom meta keys.
    def per_job = (folds_per_job as Integer) ?: 0

    SATIVAEPANG_LOOPLACE(
        SATIVAEPANG_LOOTASKS.out.taskdir
            .flatMap { meta, taskdir ->
                // 0 when the manifest is unreadable, which is what a stub run writes:
                // that falls through to one job below, same as not asking for a split.
                def manifest = taskdir.resolve('manifest.json')
                def n_folds = (manifest.size() > 0 ? new groovy.json.JsonSlurper()
                    .parseText(manifest.text).n_folds ?: 0 : 0) as Integer
                def n_shards = per_job > 0 && n_folds > per_job
                    ? (n_folds + per_job - 1).intdiv(per_job)
                    : 1
                (0..<n_shards).collect { i ->
                    [ meta + [ shard: "${i * per_job}-${Math.min((i + 1) * per_job, n_folds) - 1}", nshards: n_shards ], taskdir ]
                }
            }
    )

    def ch_placed = SATIVAEPANG_LOOPLACE.out.taskdir
        .map { meta, taskdir -> [ meta - meta.subMap('shard', 'nshards'), meta.nshards, taskdir ] }
        .branch { _meta, nshards, _taskdir ->
            one:  nshards <= 1
            many: nshards > 1
        }

    // groupKey releases each group as soon as its own shards arrive.
    SATIVAEPANG_MERGEFOLDS(
        ch_placed.many
            .map { meta, nshards, taskdir -> [ groupKey(meta, nshards), taskdir ] }
            .groupTuple()
            .map { key, taskdirs -> [ key.target, taskdirs ] }
    )

    def ch_taskdir = ch_placed.one
        .map { meta, _nshards, taskdir -> [ meta, taskdir ] }
        .mix(SATIVAEPANG_MERGEFOLDS.out.taskdir)

    // ── Phase 3: Score and report ───────────────────────────────────────────────
    //
    // refjson and the placed taskdir come from opposite ends of the chain, with no
    // other guaranteed correlation -- join explicitly rather than relying on
    // emission order.
    SATIVAEPANG_LOOSCORE(
        SATIVAEPANG_REFERENCE.out.refjson.join(ch_taskdir)
    )

    SATIVAEPANG_MISREPORT(
        SATIVAEPANG_LOOSCORE.out.mis.join(ch_taskdir)
    )

    emit:
    mislabels = SATIVAEPANG_MISREPORT.out.mislabels // [ meta, tsv ]  putative mislabels, ranked
    summary   = SATIVAEPANG_MISREPORT.out.summary   // [ meta, txt ]  run statistics
    tree      = ch_tree                            // [ meta, nwk ]  reference tree (cache for reuse)
    model     = ch_model                           // [ meta, txt ]  RAxML-NG model  (cache for reuse)
}
