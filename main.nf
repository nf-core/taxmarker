#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/taxmarker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/taxmarker
    Website: https://nf-co.re/taxmarker
    Slack  : https://nfcore.slack.com/channels/taxmarker
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TAXMARKER               } from './workflows/taxmarker'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_taxmarker_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_taxmarker_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER TYPES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Types only: the defaults live in nextflow.config
params {
    skip_clustering: Boolean
    skip_raxtax: Boolean
    skip_gapfilter: Boolean
    skip_profile_cover: Boolean
    skip_sativa: Boolean
    raxmlng_fast: Boolean
    raxmlng_seed: Integer
    version: Boolean
    plaintext_email: Boolean
    monochrome_logs: Boolean
    validate_params: Boolean
    help_full: Boolean
    show_hidden: Boolean
    min_weight: Float?
    cluster_identity: Float
    raxtax_filter_rank: Integer
    raxtax_min_confidence: Float
    min_nongap: Float
    min_profile_cover: Float
    min_profile_cover_rescue: Float
    folds_per_job: Integer?
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_TAXMARKER {

    take:
    taxonomy  // channel: taxonomy file
    sequences // channel: sequences file, aligned or not

    main:

    //
    // WORKFLOW: Run pipeline
    //
    TAXMARKER (
        taxonomy,
        sequences,
        params.seqgrep,
        params.upstream,
        params.gtdb_bac120_metadata,
        params.gtdb_ar53_metadata,
        params.skip_clustering,
        params.sequence_weights,
        params.min_weight,
        params.skip_raxtax,
        params.skip_gapfilter,
        params.skip_profile_cover,
        params.skip_sativa,
        params.taxcode,
        params.raxmlng_model,
        params.folds_per_job,
        params.export_n_per_species,
        params.hmm,
        params.hmm_name,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    multiqc_report = TAXMARKER.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.taxonomy,
        params.sequences,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_TAXMARKER (
        PIPELINE_INITIALISATION.out.taxonomy,
        PIPELINE_INITIALISATION.out.sequences
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_TAXMARKER.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
