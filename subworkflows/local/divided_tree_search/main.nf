include { TREEDIVIDE       } from '../../../modules/local/treedivide/main'
include { TREEGRAFT        } from '../../../modules/local/treegraft/main'
include { RAXMLNG_EVALUATE } from '../../../modules/local/raxmlng/evaluate/main'
include { RAXMLNG_SEARCH   } from '../../../modules/nf-core/raxmlng/search/main'

workflow DIVIDED_TREE_SEARCH {

    take:
    ch_input // channel: [ val(meta), path(alignment), path(guide_tree), val(model) ]

    main:
    TREEDIVIDE(
        ch_input.map { meta, alignment, guide_tree, _model -> [ meta, alignment, guide_tree, [] ] }
    )

    // fasta and constraints hold one file per part with the same base name; the part id
    // is the meta id so that the search names its tree <id>.raxml.bestTree, which is how
    // TREEGRAFT finds the trees again.
    def ch_parts = TREEDIVIDE.out.fasta
        .join(TREEDIVIDE.out.constraints)
        .join(ch_input.map { meta, _alignment, _guide_tree, model -> [ meta, model ] })
        .flatMap { meta, fastas, constraints, model ->
            def by_name = [ constraints ].flatten().collectEntries { c -> [ c.baseName, c ] }
            def fasta_list = [ fastas ].flatten()
            fasta_list.collect { f ->
                [ [ id: f.baseName, parent: meta.id ], f, by_name[f.baseName], model ]
            }
        }
    // ch_parts: [ part meta, fasta, constraint, model ]

    // The fasta and constraint inputs are paired by position: both come from ch_parts, so keep it that way.
    RAXMLNG_SEARCH(
        ch_parts.map { part, fasta, _constraint, model -> [ part, fasta, model ] },
        [],
        ch_parts.map { _part, _fasta, constraint, _model -> constraint },
        [],
        []
    )

    // A part whose constraint tree is fully resolved is not searched; its tree goes to the graft as it is.
    // The plan says how many trees an input has, so groupKey releases it as soon as they are all there,
    // also when none of them was searched.
    def ch_part_count = TREEDIVIDE.out.plan
        .map { meta, plan -> [ meta.id, new groovy.json.JsonSlurper().parseText(plan.text).jobs.size() ] }

    def ch_part_trees = RAXMLNG_SEARCH.out.phylogeny
        .map { part, tree -> [ part.parent, tree ] }
        .mix(TREEDIVIDE.out.resolved.flatMap { meta, trees -> [ trees ].flatten().collect { tree -> [ meta.id, tree ] } })
        .combine(ch_part_count, by: 0)
        .map { id, tree, n -> [ groupKey(id, n), tree ] }
        .groupTuple()
        .map { key, trees -> [ key.target, trees ] }

    TREEGRAFT(
        TREEDIVIDE.out.plan
            .map { meta, plan -> [ meta.id, meta, plan ] }
            .join(ch_input.map { meta, _alignment, guide_tree, _model -> [ meta.id, guide_tree ] })
            .join(ch_part_trees)
            .map { _id, meta, plan, guide_tree, trees -> [ meta, plan, trees, guide_tree ] }
    )

    RAXMLNG_EVALUATE(
        ch_input
            .map { meta, alignment, _guide_tree, model -> [ meta, alignment, model ] }
            .join(TREEGRAFT.out.tree)
            .map { meta, alignment, model, tree -> [ meta, alignment, tree, model ] }
    )

    emit:
    tree  = RAXMLNG_EVALUATE.out.phylogeny  // channel: [ val(meta), path(tree) ]   full tree, branch lengths and model optimised
    model = RAXMLNG_EVALUATE.out.best_model // channel: [ val(meta), path(model) ]  RAxML-NG model file for tree + alignment
}
