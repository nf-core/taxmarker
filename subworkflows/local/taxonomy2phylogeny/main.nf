include { TAXONOMYTREE   } from '../../../modules/local/taxonomytree/main'
include { DIVIDED_TREE_SEARCH } from '../divided_tree_search/main'
include { RAXMLNG_SEARCH      } from '../../../modules/nf-core/raxmlng/search/main'

workflow TAXONOMY2PHYLOGENY {

    take:
    ch_taxonomy_alignment // channel: [ val(meta), path(taxonomy), path(alignment), val(raxmlng_model) ]
    divided               // value:   search the parts of the taxonomy separately and graft them, instead of the whole tree at once; needs a fixed model

    main:
    TAXONOMYTREE(
        ch_taxonomy_alignment.map { meta, taxonomy, alignment, raxmlng_model -> [ meta, taxonomy ] }
    )

    // RAXMLNG_SEARCH's tree_constraint is a bare path, not a [meta, path] tuple, so it
    // can't be paired with the right alignment by a plain .map() off TAXONOMYTREE's own
    // output -- .join() first, on the shared meta, then re-project the single joined
    // channel into RAXMLNG_SEARCH's separate positional inputs. That keeps both correctly
    // paired regardless of the two processes' relative completion order under concurrency.
    def ch_search_input = ch_taxonomy_alignment
        .map { meta, taxonomy, alignment, raxmlng_model -> [ meta, alignment, raxmlng_model ] }
        .join(TAXONOMYTREE.out.guide_tree)
    // ch_search_input: [ meta, alignment, raxmlng_model, guide_tree ]

    def ch_tree
    def ch_model
    if (divided) {
        DIVIDED_TREE_SEARCH(
            ch_search_input.map { meta, alignment, raxmlng_model, guide_tree -> [ meta, alignment, guide_tree, raxmlng_model ] }
        )
        ch_tree  = DIVIDED_TREE_SEARCH.out.tree
        ch_model = DIVIDED_TREE_SEARCH.out.model
    } else {
        RAXMLNG_SEARCH(
            ch_search_input.map { meta, alignment, raxmlng_model, _guide_tree -> [ meta, alignment, raxmlng_model ] },
            [],
            ch_search_input.map { meta, alignment, raxmlng_model, guide_tree -> guide_tree },
            [],
            []
        )
        ch_tree  = RAXMLNG_SEARCH.out.phylogeny
        ch_model = RAXMLNG_SEARCH.out.best_model
    }

    emit:
    tree       = ch_tree                             // channel: [ val(meta), path(tree) ]       ML phylogeny, consistent with the taxonomy constraint
    model      = ch_model                            // channel: [ val(meta), path(model) ]      RAxML-NG model file, consistent with tree + alignment
    guide_tree = TAXONOMYTREE.out.guide_tree         // channel: [ val(meta), path(guide_tree) ] the multifurcating taxonomy-only guide tree
}
