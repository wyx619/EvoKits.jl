module EvoKits

using DataFrames
using Distributions
using LinearAlgebra
using NewickTree
using Optim
using Random
using Statistics

export CompactTree
export load_newick_tree
export to_compact_tree
export from_compact_tree
export save_newick_tree
export preprocess_tree
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
include("tools/common.jl")
include("tools/pagel_lambda.jl")
include("tools/evol_distinct.jl")
include("tools/make_ultrametric.jl")
include("tools/multi2di.jl")
include("tools/ses_pd.jl")
include("tools/ses_mpd.jl")
include("tools/ses_mntd.jl")

end
