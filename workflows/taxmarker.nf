/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_taxmarker_pipeline'
include { RESOLVETAXONOMY        } from '../modules/local/resolvetaxonomy/main'
include { GTDBWEIGHTTRANSLATE    } from '../modules/local/gtdbweighttranslate/main'
include { SEQKIT_GREP            } from '../modules/nf-core/seqkit/grep/main'
include { CHECKNAMECONSISTENCY   } from '../modules/local/checknameconsistency/main'
include { EMBOSS_SEQRET          } from '../modules/nf-core/emboss/seqret/main'
include { WEIGHTED_CLUSTERING    } from '../subworkflows/local/weighted_clustering'
include { ENSURE_ALIGNED         } from '../subworkflows/local/ensure_aligned'
include { GAPFILTER              } from '../modules/local/gapfilter/main'
include { PROFILECOVER           } from '../modules/local/profilecover/main'
include { RAXTAX_PREFILTER       } from '../subworkflows/local/raxtax_prefilter'
include { SATIVA as SWF_SATIVA   } from '../subworkflows/local/sativa'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TAXMARKER {

    take:
    ch_taxonomy        // channel: taxonomy file, or [] if not provided (derived from --sequences headers instead)
    ch_sequences       // channel: sequences file, aligned or not
    seqgrep              // value: keep only --sequences records whose header matches this pattern, or null/empty to skip filtering
    upstream             // value: name of an upstream data source ('GTDB'), or null/empty for none
    gtdb_bac120_metadata // value: path to GTDB bac120 metadata, or null/empty if not needed
    gtdb_ar53_metadata   // value: path to GTDB ar53 metadata, or null/empty if not needed
    skip_clustering      // value: skip WEIGHTED_CLUSTERING entirely?
    sequence_weights   // value:   path to a two-column weight table, or null/empty if not supplied
    min_weight         // value:   absolute weight cutoff, or null/empty to skip it
    skip_raxtax        // value:   skip the raxtax prefilter?
    skip_gapfilter     // value:   skip the gap filter (already-aligned input)?
    skip_profile_cover // value:   skip the profile-coverage filter (hmmalign-derived input)?
    skip_sativa        // value:   skip the phylogenetic placement subworkflow entirely?
    taxcode            // value:   taxonomic code for sativa-epang (bac/bot/zoo/vir)
    placement_jobs     // value:   jobs the leave-one-out placement is spread over
    hmm                // value:   path to an HMM profile database, or null/empty if not needed
    hmm_name           // value:   name of a specific profile within hmm, or null/empty
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'taxmarker_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: SEQKIT_GREP (optional, --seqgrep to enable)
    //
    // Keep only matching --sequences records before anything else runs, e.g. to
    // restrict a combined multi-taxon reference set to the records wanted.
    //
    def ch_sequences_for_resolve = ch_sequences
    if (seqgrep) {
        SEQKIT_GREP(ch_sequences.map { [ [ id: 'user-alignment' ], it ] }, [], '')
        ch_sequences_for_resolve = SEQKIT_GREP.out.filter.map { _meta, fasta -> fasta }
    }

    //
    // MODULE: RESOLVETAXONOMY
    //
    // Resolve taxonomy from an explicit --taxonomy file if given; otherwise derive it
    // from --sequences record headers instead (GTDB-style: >id taxonomy;string), no
    // separate mode-switch param needed. Headers are always stripped down to a bare
    // id either way. If both a file and embedded header text are present, the file
    // wins -- warned about, not silently ignored.
    //
    RESOLVETAXONOMY(
        // '' (no taxonomy given, see PIPELINE_INITIALISATION) becomes [] here, right
        // at the input tuple RESOLVETAXONOMY itself receives -- the literal empty
        // list Nextflow recognises as "optional path input, absent" when it's one
        // element of a freshly-built tuple, as opposed to a channel item in its own
        // right (which .combine() would silently flatten away).
        ch_taxonomy.combine(ch_sequences_for_resolve).map { tax, seq -> [ [ id: 'user-alignment' ], tax ?: [], seq ] }
    )
    RESOLVETAXONOMY.out.warnings.subscribe { _meta, warnings_file ->
        def text = warnings_file.text.trim()
        if (text) {
            log.warn(text)
        }
    }
    def ch_taxonomy_resolved  = RESOLVETAXONOMY.out.taxonomy.map  { _meta, tax -> tax }
    def ch_sequences_resolved = RESOLVETAXONOMY.out.sequences.map { _meta, seq -> seq }

    //
    // MODULE: Validate that taxonomy and sequences name the same records, and
    // rewrite characters that are difficult for downstream tools (e.g. parens) in
    // both. Runs first, as a process (not inline Nextflow code) so a large input
    // doesn't inflate the head job's memory/CPU footprint.
    //
    CHECKNAMECONSISTENCY(
        ch_taxonomy_resolved.combine(ch_sequences_resolved).map { tax, seq -> [ [ id: 'user-alignment' ], tax, seq ] }
    )
    def ch_taxonomy_checked  = CHECKNAMECONSISTENCY.out.checked.map { _meta, tax, _seq -> tax }
    def ch_sequences_checked = CHECKNAMECONSISTENCY.out.checked.map { _meta, _tax, seq -> seq }

    //
    // MODULE: Normalise the sequences to FASTA once, here, rather than separately
    // inside RAXTAX_PREFILTER and SWF_SATIVA (previously duplicated). Also gives
    // ENSURE_ALIGNED a single canonical format to inspect for the unaligned-input
    // support below.
    //
    EMBOSS_SEQRET(ch_sequences_checked.map { [ [ id: 'user-alignment' ], it ] }, 'fasta')
    def ch_sequences_fasta = EMBOSS_SEQRET.out.outseq.map { _meta, seq -> seq }

    //
    // MODULE: GTDBWEIGHTTRANSLATE (optional, --upstream 'GTDB' to enable)
    //
    // Translates GTDB genome metadata into the generic --sequence_weights table
    // WEIGHTED_CLUSTERING consumes below -- see nf-core/taxmarker#15. Runs on the
    // pristine, pre-RESOLVETAXONOMY sequences: it needs both the GTDB accession
    // (before CHECKNAMECONSISTENCY sanitises '~' away) and each record's own header
    // metadata (e.g. [contig_len=...]), neither of which survive past this point.
    //
    // Every branch below wraps its value in a 1-element list -- see WEIGHTED_CLUSTERING's
    // own comment on ch_sequence_weights for why: .combine() would otherwise flatten a
    // bare [] payload's own (zero) elements in, instead of one empty-list slot.
    def ch_sequence_weights
    if (upstream == 'GTDB') {
        if (!gtdb_bac120_metadata && !gtdb_ar53_metadata) {
            error("--upstream 'GTDB' needs --gtdb_bac120_metadata and/or --gtdb_ar53_metadata.")
        }
        if (sequence_weights) {
            // Same precedent as --taxonomy winning over embedded taxonomy text.
            log.warn("Both --upstream 'GTDB' and --sequence_weights were given; the explicit --sequence_weights file wins.")
            ch_sequence_weights = channel.value([ file(sequence_weights, checkIfExists: true) ])
        } else {
            GTDBWEIGHTTRANSLATE(
                ch_sequences_for_resolve.map { [ [ id: 'user-alignment' ], it ] },
                gtdb_bac120_metadata ? file(gtdb_bac120_metadata, checkIfExists: true) : [],
                gtdb_ar53_metadata   ? file(gtdb_ar53_metadata, checkIfExists: true)   : []
            )
            ch_sequence_weights = GTDBWEIGHTTRANSLATE.out.weights.map { _meta, weights -> [ weights ] }
        }
    } else {
        ch_sequence_weights = channel.value([ sequence_weights ? file(sequence_weights, checkIfExists: true) : [] ])
    }

    //
    // SUBWORKFLOW: WEIGHTED_CLUSTERING (optional, skip_clustering to disable)
    //
    // Reduce the input set to one representative per (cluster, taxon) pair before
    // raxtax, alignment and placement all see it -- see nf-core/taxmarker#15.
    // --cluster_identity is read directly from conf/modules.config's ext.args, like
    // this pipeline's other tuning-only knobs.
    //
    def ch_taxonomy_clustered
    def ch_sequences_clustered
    def run_clustering = !skip_clustering.toString().toBoolean()
    if (run_clustering) {
        WEIGHTED_CLUSTERING(ch_taxonomy_checked, ch_sequences_fasta, ch_sequence_weights, min_weight)
        ch_taxonomy_clustered  = WEIGHTED_CLUSTERING.out.taxonomy
        ch_sequences_clustered = WEIGHTED_CLUSTERING.out.sequences
    } else {
        ch_taxonomy_clustered  = ch_taxonomy_checked
        ch_sequences_clustered = ch_sequences_fasta
    }

    //
    // SUBWORKFLOW: RAXTAX_PREFILTER (optional, skip_raxtax to disable)
    //
    // Fast raxtax self-classification prefilter ahead of the expensive alignment and
    // EPA-ng-based placement below. Runs on unaligned sequences (raxtax classifies
    // plain sequences, never needs an alignment) so sequences it flags skip alignment
    // entirely instead of only skipping placement. They never reach SWF_SATIVA --
    // they're reported directly via ch_raxtax_mislabels instead.
    //
    def ch_taxonomy_for_alignment
    def ch_sequences_for_alignment
    def ch_raxtax_mislabels
    // Coerce explicitly: a CLI-supplied `--skip_raxtax false` arrives as the *string*
    // "false", and Groovy's `!"false"` is false (any non-empty string is truthy) --
    // .toBoolean() parses both real Booleans and "true"/"false" strings correctly.
    // Confirmed empirically that nf-schema's cli_typecast (enabled just above, in
    // PIPELINE_INITIALISATION) validates the string against the boolean schema type but
    // does not itself replace params.skip_raxtax with a real Boolean, so this is still
    // needed even with cli_typecast on.
    def run_raxtax = !skip_raxtax.toString().toBoolean()
    if (run_raxtax) {
        RAXTAX_PREFILTER(ch_taxonomy_clustered, ch_sequences_clustered)
        ch_taxonomy_for_alignment  = RAXTAX_PREFILTER.out.taxonomy
        ch_sequences_for_alignment = RAXTAX_PREFILTER.out.sequences
        ch_raxtax_mislabels        = RAXTAX_PREFILTER.out.mislabels
    } else {
        ch_taxonomy_for_alignment  = ch_taxonomy_clustered
        ch_sequences_for_alignment = ch_sequences_clustered
        ch_raxtax_mislabels        = channel.empty()
    }

    //
    // SUBWORKFLOW: ENSURE_ALIGNED
    //
    // Transparently accepts unaligned input too, with no separate mode-switch param:
    // already-aligned content passes straight through; unaligned content is aligned
    // via hmmalign against the hmm/hmm_name profile before continuing. Only past this
    // point is the data actually guaranteed to be an alignment.
    //
    ENSURE_ALIGNED(ch_sequences_for_alignment, hmm, hmm_name)

    //
    // MODULE: GAPFILTER + PROFILECOVER (each optional, own skip flag)
    //
    // Both drop sequences too short/incomplete to place reliably, reporting them
    // separately rather than silently discarding them -- but applied to ENSURE_ALIGNED's
    // two branches independently, since they have structurally different non-gap
    // distributions (see ensure_aligned/main.nf): GAPFILTER (params.min_nongap) for
    // already-aligned input, PROFILECOVER (params.min_profile_cover) for hmmalign-
    // derived input, already masked down to the HMM's match-state columns.
    //
    def ch_taxonomy_gapfiltered
    def ch_alignment_gapfiltered
    // Coerce explicitly: a CLI-supplied `--skip_gapfilter false` arrives as the
    // *string* "false" -- see the analogous skip_raxtax coercion above for why
    // .toString().toBoolean() is needed even with nf-schema's cli_typecast enabled.
    def run_gapfilter = !skip_gapfilter.toString().toBoolean()
    if (run_gapfilter) {
        GAPFILTER(
            ch_taxonomy_for_alignment.combine(ENSURE_ALIGNED.out.alignment_passthrough).map { tax, aln -> [ [ id: 'user-alignment' ], tax, aln ] }
        )
        ch_taxonomy_gapfiltered  = GAPFILTER.out.taxonomy.map { _meta, tax -> tax }
        ch_alignment_gapfiltered = GAPFILTER.out.alignment.map { _meta, aln -> aln }
    } else {
        // ch_taxonomy_for_alignment alone would still emit its one item even when
        // alignment_passthrough is empty (e.g. unaligned input took the hmm branch
        // instead), desyncing this pair's cardinality -- 1 taxonomy item vs 0
        // alignment items -- which then corrupts the .mix() below (the stray taxonomy
        // item pairs with the *other* branch's real alignment downstream). Gate it by
        // the same alignment channel instead, so it collapses to 0 items exactly when
        // alignment_passthrough does.
        ch_alignment_gapfiltered = ENSURE_ALIGNED.out.alignment_passthrough
        ch_taxonomy_gapfiltered  = ch_taxonomy_for_alignment.combine(ch_alignment_gapfiltered).map { tax, _aln -> tax }
    }

    def ch_taxonomy_covfiltered
    def ch_alignment_covfiltered
    def run_profile_cover = !skip_profile_cover.toString().toBoolean()
    if (run_profile_cover) {
        PROFILECOVER(
            ch_taxonomy_for_alignment.combine(ENSURE_ALIGNED.out.alignment_from_hmm).map { tax, aln -> [ [ id: 'user-alignment' ], tax, aln ] }
        )
        ch_taxonomy_covfiltered  = PROFILECOVER.out.taxonomy.map { _meta, tax -> tax }
        ch_alignment_covfiltered = PROFILECOVER.out.alignment.map { _meta, aln -> aln }
    } else {
        // See the analogous gapfilter comment above -- same cardinality-gating fix.
        ch_alignment_covfiltered = ENSURE_ALIGNED.out.alignment_from_hmm
        ch_taxonomy_covfiltered  = ch_taxonomy_for_alignment.combine(ch_alignment_covfiltered).map { tax, _aln -> tax }
    }

    // Exactly one of ENSURE_ALIGNED's two branches ever has content for a given run
    // (CHECKALIGNED classifies the whole input as aligned-or-not, never a mix), so
    // .mix() here just recombines whichever branch actually ran with the other's
    // always-empty channel.
    def ch_taxonomy_for_sativa  = ch_taxonomy_gapfiltered.mix(ch_taxonomy_covfiltered)
    def ch_alignment_for_sativa = ch_alignment_gapfiltered.mix(ch_alignment_covfiltered)

    //
    // SUBWORKFLOW: SATIVA (optional, skip_sativa to disable)
    //
    // This implements all the logic in the workflow: builds the reference tree,
    // leave-one-out places every sequence back into it, and scores each one.
    // Skipping it turns the rest of the pipeline into a general-purpose taxonomy-
    // resolution/alignment/prefilter QC tool -- e.g. to get a cleaned, filtered
    // alignment+taxonomy pair (already published by whichever upstream module
    // produced it last) or raxtax-only mislabels, without the much more expensive
    // EPA-ng-based placement step.
    //
    // The later two SWF_SATIVA params are meant to pass a reference tree and a
    // model file respectively. Not implemented yet.
    //
    def ch_sativa_mislabels
    def run_sativa = !skip_sativa.toString().toBoolean()
    if (run_sativa) {
        SWF_SATIVA(ch_taxonomy_for_sativa, ch_alignment_for_sativa, taxcode, placement_jobs, [], [])
        ch_sativa_mislabels = SWF_SATIVA.out.mislabels
    } else {
        ch_sativa_mislabels = channel.empty()
    }

    //
    // Merge raxtax-flagged mislabels (skipped placement entirely) with SATIVASCORE's own
    // into one final report. Both share the same TSV schema (method column distinguishes
    // detection source), so collectFile with keepHeader can concatenate them directly --
    // no bridging process needed just to reshape/combine two files.
    //
    ch_raxtax_mislabels
        .mix(ch_sativa_mislabels)
        .map { _meta, tsv -> tsv }
        .collectFile(name: 'user-alignment.mislabels.tsv', storeDir: "${outdir}/mislabels", keepHeader: true, skip: 1)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'taxmarker'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
