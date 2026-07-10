module EvoKits

using MKL
using CSV
using DataFrames
using Distributions
using LinearAlgebra
using NewickTree
using NLopt
using Optim
using Random
using Statistics

export CompactTree
export load_newick_tree
export to_compact_tree
export from_compact_tree
export save_newick_tree
export PhyloRef
export build_phyloref
export phyloref_node_table
export phyloref_edge_table
export phylo_node_anchor
export phylo_branch_anchor
export phylo_edge_signature
export phylo_edge_signatures
export phylo_edges_from_signatures
export ape_node_id
export compact_node_id
export ape_edge_id
export compact_edge_id
export ape_edge_matrix
export attach_ape_node_map
export attach_ape_edge_map
export align_traits_to_tree
export compute_aic
export compute_aicc
export compute_bic
export pagel_lambda
export summary_pagel_lambda
export EvolutionaryDistinctivenessResult
export evol_distinct
export evol_distinct_dataframe
export make_ultrametric
export multi2di
export keep_tip
export drop_tip
export pd
export pd_dataframe
export ses_pd
export ses_pd_dataframe
export ses_mpd
export ses_mpd_dataframe
export ses_mntd
export ses_mntd_dataframe

include("types.jl")
include("io.jl")
include("checks.jl")
include("criteria.jl")
include("trait_alignment.jl")
include("phyloref.jl")
include("tools/common.jl")
include("tools/pagel_lambda.jl")
include("tools/evol_distinct.jl")
include("tools/make_ultrametric.jl")
include("tools/multi2di.jl")
include("tree_pruning.jl")
include("tools/ses_pd.jl")
include("tools/ses_mpd.jl")
include("tools/ses_mntd.jl")

end
