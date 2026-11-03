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
using SparseArrays
using Statistics
using StableRNGs

export CompactTree
export read_tree
export write_tree
export mrca
export distance
export pairwise_distance
export distances_from_root
export distances_to_tip
export descendant_tips
export clade_members
export subtree_at_node
export CladeGroupingResult
export get_clades
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
export timeslice
export multi2di
export keep_tip
export drop_tip
export pd
export pd_dataframe
export WeightedEndemismResult
export weighted_endemism
export PhyloEndemismResult
export phylo_endemism
export RelativePDResult
export relative_pd
export PhylobetaResult
export phylobeta
export BetaDissResult
export beta_diss
export CommunityContext
export prepare_community
export CommunityData
export align_community
export CommunityTraitResult
export map_trait
export aggregate_species
export LineageCommunityResult
export collapse_range
export FunctionalGroupingResult
export functional_groups
export collapse_traits
export FunctionalBetaResult
export functional_beta
export FunctionalGroupSelectionResult
export optimal_functional_groups
export IndicatorSpeciesResult
export indicators
export ValueClassificationResult
export classify_values
export PhyloBuilderResult
export phylobuilder
export EDGEResult
export edge
export MetricSelectionResult
export hotspots
export coldspots
export MeanDistResult
export mean_dist
export PhylobetaSESResult
export phylobeta_ses
export PhyloregionResult
export phyloregion
export PhyloregionContext
export prepare_phyloregion
export cutree
export LinkageSelectionResult
export select_linkage
export OptimalPhyloregionResult
export optimal_phyloregion
export NMDSResult
export nmds
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
export nri_nti_table
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

include("base/labels.jl")
include("base/criteria.jl")
include("base/random.jl")
include("base/tree_types.jl")
include("base/tree_io.jl")

include("tree/queries/validation.jl")
include("tree/queries/query.jl")
include("tree/queries/distances.jl")
include("tree/queries/clades.jl")
include("tree/queries/clade_groups.jl")
include("tree/operations/pruning.jl")
include("tree/mapping/mapping.jl")
include("tree/operations/transforms.jl")
include("tree/operations/multi2di.jl")
include("tree/operations/timeslice.jl")
include("tree/metrics/metrics.jl")
include("tree/metrics/shape.jl")
include("tree/operations/phylobuilder.jl")

include("traits/alignment.jl")
include("traits/pagel_lambda.jl")
include("traits/phylo_d.jl")
include("traits/signal_k.jl")
include("traits/pic.jl")

include("community/core/input.jl")
include("community/core/branches.jl")
include("community/null_models/null_models.jl")
include("community/diversity/pd.jl")
include("community/diversity/endemism.jl")
include("community/diversity/mpd.jl")
include("community/diversity/mntd.jl")
include("community/diversity/beta.jl")
include("community/diversity/taxonomic_beta.jl")
include("community/core/context.jl")
include("community/diversity/diversity.jl")
include("community/core/data_metrics.jl")
include("community/aggregation/aggregation.jl")
include("community/aggregation/lineage_aggregation.jl")
include("community/aggregation/functional_beta.jl")
include("community/aggregation/functional_selection.jl")
include("community/regionalization/indicator_species.jl")
include("community/conservation/conservation.jl")
include("community/conservation/value_classification.jl")
include("community/regionalization/regionalization.jl")
include("community/regionalization/regionalization_context.jl")
include("community/regionalization/regionalization_selection.jl")
include("community/regionalization/mean_dist.jl")
include("community/regionalization/nmds.jl")
include("community/null_models/nri_nti.jl")
include("community/null_models/ses.jl")
include("community/null_models/beta_ses.jl")
include("community/null_models/nri_nti_table.jl")
include("community/conservation/edge2.jl")

end
