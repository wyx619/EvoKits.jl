module EvoKits

# MKL must load first to select the BLAS backend.
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

# Shared numerical helpers
export compute_aic, compute_aicc, compute_bic

# Tree representation and I/O
export CompactTree, read_tree, write_tree

# Tree queries, distances and clade grouping
export mrca, descendant_tips, clade_members, subtree_at_node
export distance, pairwise_distance, distances_from_root, distances_to_tip, node_depths
export CladeGroupingResult, get_clades

# Tree transformations and completion
export keep_tip, drop_tip, make_ultrametric, multi2di, timeslice
export PhyloBuilderResult, phylobuilder

# Node and edge identities, including R/ape numbering
export PhyloRef, build_phyloref, phyloref_node_table, phyloref_edge_table
export phylo_node_anchor, phylo_branch_anchor
export phylo_edge_signature, phylo_edge_signatures, phylo_edges_from_signatures
export ape_node_id, compact_node_id, ape_edge_id, compact_edge_id
export ape_edge_matrix
export attach_ape_node_map, attach_ape_edge_map

# Tree metrics
export EvolutionaryDistinctivenessResult, evol_distinct, evol_distinct_dataframe
export tree_imbalance, gamma_statistic
export LineageThroughTimeResult, lineages_through_time

# Trait alignment and phylogenetic signal
export align_traits_to_tree
export pagel_lambda, summary_pagel_lambda
export PhyloDResult, phylo_d, PhyloDSubsetResult, phylo_d_subset
export PhyloSignalResult, phylogenetic_signal_k
export PICResult, phylogenetic_contrasts

# Community input and reusable contexts
export CommunityData
export align_community
export CommunityContext, prepare_community

# Community diversity and endemism
export pd, pd_dataframe
export WeightedEndemismResult, weighted_endemism
export PhyloEndemismResult, phylo_endemism
export RelativePDResult, relative_pd
export RaoDResult, raoD
export PSVResult, psv, psr, pse, psc, psd
export ExpectedPDResult, expected_pd, pd_variance
export edge_abundance_distribution, ead, pd_calc
export PDBootstrapResult, pd_bootstrap

# Pairwise and multiple-site beta diversity
export PhylobetaResult, phylobeta
export BetaDissResult, beta_diss
export BetaMultiResult, beta_multi
export PhylobetaMultiResult, phylobeta_multi
export comdist, comdistnt, unifrac, phylosor

# Community aggregation and functional beta diversity
export CommunityTraitResult, map_trait
export aggregate_species
export LineageCommunityResult, collapse_range
export FunctionalGroupingResult, functional_groups, collapse_traits
export FunctionalBetaResult, functional_beta
export FunctionalGroupSelectionResult, optimal_functional_groups

# Regionalization, mixed membership and indicator species
export PhyloregionResult, phyloregion
export PhyloregionContext, prepare_phyloregion
export cutree
export LinkageSelectionResult, select_linkage
export OptimalPhyloregionResult, optimal_phyloregion
export GoMResult, fitgom
export IndicatorSpeciesResult, indicators

# Distance aggregation and ordination
export MeanDistResult, mean_dist
export NMDSResult, nmds
export NMDSAlignmentResult, align_nmds
export NMDSColorScale, NMDSColorResult, nmds_colors

# Null models and community structure
export stable_grid_seed
export ses_pd, ses_pd_dataframe
export ses_mpd, ses_mpd_dataframe
export ses_mntd, ses_mntd_dataframe
export PhylobetaSESResult, phylobeta_ses
export NriNtiResult, nri_nti, nri_nti_table

# Conservation scores and value classification
export EDGEResult, edge
export MetricSelectionResult, hotspots, coldspots
export ValueClassificationResult, classify_values
export Edge2Result, edge2_mod, score_adjust
export PextDistribution, calc_pext_distribution, load_pext_distribution_csv
export sample_pext, reorder_pext_to_tree

# Source loading follows dependencies; directory groups may span several phases.

# Shared helpers and tree representation
include("base/labels.jl")
include("base/criteria.jl")
include("base/random.jl")
include("base/tree_types.jl")
include("base/tree_io.jl")

# Tree queries, operations, mapping and metrics
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

# Trait alignment and phylogenetic signal
include("traits/alignment.jl")
include("traits/pagel_lambda.jl")
include("traits/phylo_d.jl")
include("traits/signal_k.jl")
include("traits/pic.jl")

# Community input, branch traversal and randomization primitives
include("community/core/input.jl")
include("community/core/branches.jl")
include("community/null_models/null_models.jl")

# Community diversity kernels
include("community/diversity/pd.jl")
include("community/diversity/endemism.jl")
include("community/diversity/mpd.jl")
include("community/diversity/mntd.jl")
include("community/diversity/beta.jl")
include("community/diversity/taxonomic_beta.jl")
include("community/diversity/beta_multi.jl")
include("community/diversity/phylobeta_multi.jl")

# Reusable community context and metric adapters
include("community/core/context.jl")
include("community/diversity/diversity.jl")
include("community/core/data_metrics.jl")

# Species, lineage and functional aggregation
include("community/aggregation/aggregation.jl")
include("community/aggregation/lineage_aggregation.jl")
include("community/aggregation/functional_beta.jl")
include("community/aggregation/functional_selection.jl")

# Indicator species and fixed-k mixture models
include("community/regionalization/indicator_species.jl")
include("community/regionalization/gom_fit.jl")
include("community/regionalization/gom.jl")

# Conservation scores and value classification
include("community/conservation/conservation.jl")
include("community/conservation/value_classification.jl")

# Regionalization and ordination
include("community/regionalization/regionalization.jl")
include("community/regionalization/regionalization_context.jl")
include("community/regionalization/regionalization_selection.jl")
include("community/regionalization/mean_dist.jl")
include("community/regionalization/nmds.jl")
include("community/regionalization/nmds_alignment.jl")
include("community/regionalization/nmds_colors.jl")

# Null-model statistics and community structure
include("community/null_models/nri_nti.jl")
include("community/null_models/ses.jl")
include("community/null_models/beta_ses.jl")
include("community/null_models/nri_nti_table.jl")

# EDGE2 conservation assessment
include("community/conservation/edge2.jl")

end
