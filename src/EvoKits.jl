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
using StableRNGs

export CompactTree
export load_newick_tree
export to_compact_tree
export from_compact_tree
export save_newick_tree
export mrca
export distance
export pairwise_distance
export distances_from_root
export distances_to_tip
export descendant_tips
export clade_members
export subtree_at_node
export node_depths
export tree_imbalance
export gamma_statistic
export LineageThroughTimeResult
export lineages_through_time
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
export PhyloDResult
export phylo_d
export PhyloDSubsetResult
export phylo_d_subset
export Edge2Result
export PextDistribution
export edge2_mod
export calc_pext_distribution
export load_pext_distribution_csv
export sample_pext
export reorder_pext_to_tree
export score_adjust
export NriNtiResult
export nri_nti
export stable_grid_seed
export comdist
export comdistnt
export unifrac
export phylosor
export RaoDResult
export raoD
export PSVResult
export psv
export psr
export pse
export psc
export psd
export ExpectedPDResult
export expected_pd
export pd_variance
export edge_abundance_distribution
export ead
export pd_calc
export PDBootstrapResult
export pd_bootstrap
export PhyloSignalResult
export phylogenetic_signal_k
export PICResult
export phylogenetic_contrasts

include("base/criteria.jl")
include("base/random.jl")
include("base/tree_types.jl")
include("base/tree_io.jl")

include("tree/validation.jl")
include("tree/query.jl")
include("tree/distances.jl")
include("tree/clades.jl")
include("tree/pruning.jl")
include("tree/mapping.jl")
include("tree/transforms.jl")
include("tree/multi2di.jl")
include("tree/metrics.jl")
include("tree/shape.jl")

include("traits/alignment.jl")
include("traits/pagel_lambda.jl")
include("traits/phylo_d.jl")
include("traits/signal_k.jl")
include("traits/pic.jl")

include("community/input.jl")
include("community/null_models.jl")
include("community/pd.jl")
include("community/mpd.jl")
include("community/mntd.jl")
include("community/beta.jl")
include("community/diversity.jl")
include("community/nri_nti.jl")
include("community/edge2.jl")

end
