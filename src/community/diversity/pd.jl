Base.@kwdef struct PDResult
    sample_labels::Vector{String} = String[]
    pd::Vector{Float64} = Float64[]
    richness::Vector{Int} = Int[]
    include_root::Bool = true
end

Base.@kwdef struct SESPDResult
    sample_labels::Vector{String} = String[]
    ntaxa::Vector{Int} = Int[]
    pd_obs::Vector{Float64} = Float64[]
    pd_rand_mean::Vector{Float64} = Float64[]
    pd_rand_sd::Vector{Float64} = Float64[]
    pd_obs_rank::Vector{Float64} = Float64[]
    pd_obs_z::Vector{Float64} = Float64[]
    pd_obs_p::Vector{Float64} = Float64[]
    runs::Int = 0
    null_model::Symbol = :unknown
    include_root::Bool = true
end


function pd(samp::AbstractMatrix{<:Real}, tree::CompactTree; include_root::Bool = true, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing)
    _validate_community_branch_lengths(tree)
    labels = _tools_sample_labels(size(samp, 1), sample_labels)
    rows = _community_tip_rows(samp, tree, tip_labels)
    richness = length.(rows)
    roots = include_root ? nothing : _community_roots(rows, tree)
    pd_values = zeros(Float64, length(rows))
    _walk_community_branches(rows, tree; roots = roots) do site, edge
        pd_values[site] += tree.edge_length[edge]
    end
    if !include_root
        pd_values[richness .== 1] .= NaN
    end
    return PDResult(sample_labels = labels, pd = pd_values, richness = richness, include_root = include_root)
end

"""Relative PD and its normalized original/equal-positive-branch components."""
struct RelativePDResult
    sample_labels::Vector{String}
    rpd::Vector{Float64}
    pd_normalized::Vector{Float64}
    pd_equal_normalized::Vector{Float64}
    richness::Vector{Int}
    total_branch_length::Float64
    positive_branch_count::Int
end

"""
    relative_pd(community, tree; tip_labels=nothing, sample_labels=nothing)

Divide PD on the original tree by PD on a tree with the same topology and
equal positive branch lengths. Normalize each tree's total branch length to
one; zero branches remain zero. Includes root paths and never modifies tree.
Zero denominator gives NaN, including empty sites and all-zero-length trees.
"""
function relative_pd(community::AbstractMatrix{<:Real}, tree::CompactTree;
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    _validate_community_branch_lengths(tree)
    labels = _tools_sample_labels(size(community, 1), sample_labels)
    rows = _community_tip_rows(community, tree, tip_labels)
    total = sum(tree.edge_length)
    isfinite(total) || throw(ArgumentError("total branch length must be finite"))
    positive = count(>(0), tree.edge_length)
    original, equal = zeros(Float64, length(rows)), zeros(Float64, length(rows))
    _walk_community_branches(rows, tree) do site, edge
        original[site] += tree.edge_length[edge]
        tree.edge_length[edge] > 0 && (equal[site] += 1)
    end
    if total > 0
        original ./= total
        equal ./= positive
    else
        fill!(original, NaN)
        fill!(equal, NaN)
    end
    rpd = [equal[i] > 0 ? original[i] / equal[i] : NaN for i in eachindex(equal)]
    return RelativePDResult(labels, rpd, original, equal, length.(rows), total, positive)
end

function pd_dataframe(samp::AbstractMatrix{<:Real}, tree::CompactTree; kwargs...)
    result = pd(samp, tree; kwargs...)
    return DataFrame(sample = result.sample_labels, PD = result.pd, SR = result.richness)
end

"""
    ses_pd(community, tree; null_model=:taxa_labels, runs=999, iterations=1000,
           include_root=true, tip_labels=nothing, sample_labels=nothing,
           rng=Random.default_rng(), max_attempts_per_swap=10_000)

Standardized effect size of rooted PD. Observed and randomized communities use
sparse species sets; no dense site-by-species buffers or all-tip root paths are
stored. Null means and sample standard deviations use finite draws and Welford
updates. Sites with no finite draws have NaN mean/SD/z; zero null SD gives NaN z.
Rank p-values retain (1 + count_below + half_ties) / (runs + 1).
"""
function ses_pd(samp::AbstractMatrix{<:Real}, tree::CompactTree; null_model::Symbol = :taxa_labels, runs::Integer = 999, iterations::Integer = 1000, include_root::Bool = true, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, rng::AbstractRNG = Random.default_rng(), max_attempts_per_swap::Integer = 10_000)
    result = _ses_summary(:pd, samp, tree; null_model, runs, iterations, include_root,
        tip_labels, sample_labels, rng, max_attempts_per_swap)
    return SESPDResult(sample_labels = result.labels, ntaxa = result.richness,
        pd_obs = result.obs, pd_rand_mean = result.means, pd_rand_sd = result.sds,
        pd_obs_rank = result.ranks, pd_obs_z = result.z, pd_obs_p = result.p,
        runs = Int(runs), null_model = result.model, include_root = include_root)
end

function ses_pd_dataframe(samp::AbstractMatrix{<:Real}, tree::CompactTree; kwargs...)
    result = ses_pd(samp, tree; kwargs...)
    return DataFrame(sample = result.sample_labels, ntaxa = result.ntaxa, pd_obs = result.pd_obs, pd_rand_mean = result.pd_rand_mean, pd_rand_sd = result.pd_rand_sd, pd_obs_rank = result.pd_obs_rank, pd_obs_z = result.pd_obs_z, pd_obs_p = result.pd_obs_p, runs = fill(result.runs, length(result.pd_obs)))
end
