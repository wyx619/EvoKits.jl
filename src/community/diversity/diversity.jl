Base.@kwdef struct PDBootstrapResult
    observed::Vector{Float64} = Float64[]
    replicates::Matrix{Float64} = zeros(Float64, 0, 0)
    richness::Vector{Int} = Int[]
end

function _normalize_pd_variant(variant)
    value = Symbol(lowercase(String(variant)))
    value in (:tbl, :mst, :ueh, :sbl, :tip) || throw(ArgumentError("variant must be TBL, MST, UEH, SBL, or TIP"))
    return value
end

"""
    pd_calc(community, tree; variant=:tbl, tip_labels=nothing)
    pd_calc(tree, tips; variant=:tbl)

Branch-length sums for TBL (occupied root paths), MST (connections below the
observed MRCA), UEH (one present descendant), SBL (at least two), or TIP
(original terminal edges). Presence records retain dense/CSC input semantics.
Empty sites and singleton MST give zero. No dense community or site-by-edge
matrix is built. Sparse sites use connecting trees and stable path sums;
dense sites use sequential descendant counts. Inputs are unchanged.
"""
function pd_calc(community::AbstractMatrix{<:Real}, tree::CompactTree;
    variant = :tbl, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    variant = _normalize_pd_variant(variant)
    _validate_community_branch_lengths(tree)
    rows = _community_tip_rows(community, tree, tip_labels)
    workspace = _PDWorkspace(tree)
    values = [_pd_site!(workspace, tree, tips, variant) for tips in rows]
    return (; pd = values, richness = length.(rows), variant)
end

function pd_calc(tree::CompactTree, tips; variant = :tbl)
    variant = _normalize_pd_variant(variant)
    _validate_community_branch_lengths(tree)
    positions = _pd_tip_positions(tree, tips)
    return _pd_site!(_PDWorkspace(tree), tree, positions, variant)
end

"""
    pd_bootstrap(community, tree; reps=999, tip_weights=nothing,
                 include_root=true, tip_labels=nothing, rng=Random.default_rng())
    pd_bootstrap(tree, ntips; reps=1000, method=:TBL, tip_weights=nothing,
                 rng=Random.default_rng())

Weighted sampling without replacement at the observed/fixed richness using
exponential keys. Weights follow tree-tip order. Reuses one tree/path index,
one site workspace and a heap of min(r,p-r) keys, avoiding dense randomized
communities and repeated full sorts. Ordinary finite-sum weights preserve
the preceding implementation's selected subsets and RNG consumption, with
ties resolved by original tip position. Zero-probability tips have key -Inf;
if richness exceeds positive support, index-ordered -Inf ties fill the draw,
preserving the existing convention. Extremely overflowing weight sums are
scaled before normalization. All weights must fit finite Float64.

Community output remains PDBootstrapResult with observed PD, a sites-by-reps
replicate matrix and richness. Empty sites give zero; without the root,
singleton sites give NaN. Tree-only output is a reps-length vector for the
selected PD variant. Repeated seeds reproduce draws; inputs are not mutated.
"""
function pd_bootstrap(community::AbstractMatrix{<:Real}, tree::CompactTree;
    reps::Integer = 999, tip_weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
    include_root::Bool = true, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    rng::AbstractRNG = Random.default_rng(),
)
    reps = _pd_integer(reps, "reps", 1)
    _validate_community_branch_lengths(tree)
    rows = _community_tip_rows(community, tree, tip_labels)
    probabilities = _pd_probabilities(tree, tip_weights)
    workspace, sampler = _PDWorkspace(tree), _PDSamplingWorkspace()
    richness = length.(rows)
    function evaluate(tips)
        !include_root && length(tips) == 1 && return NaN
        return _pd_site!(workspace, tree, tips, include_root ? :tbl : :mst)
    end
    observed = [evaluate(tips) for tips in rows]
    draws = Matrix{Float64}(undef, length(rows), reps)
    for rep in 1:reps, site in eachindex(rows)
        draws[site, rep] = richness[site] == 0 ? 0.0 :
            evaluate(_pd_sample!(sampler, probabilities, richness[site], rng))
    end
    return PDBootstrapResult(; observed, replicates = draws, richness)
end

function pd_bootstrap(tree::CompactTree, ntips::Integer;
    reps::Integer = 1000, method = :TBL,
    tip_weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
    rng::AbstractRNG = Random.default_rng(),
)
    ntips = _pd_integer(ntips, "ntips", 1)
    ntips <= tree.ntips || throw(ArgumentError("ntips must be between 1 and tree.ntips"))
    reps = _pd_integer(reps, "reps", 1)
    variant = _normalize_pd_variant(method)
    _validate_community_branch_lengths(tree)
    probabilities = _pd_probabilities(tree, tip_weights)
    workspace, sampler = _PDWorkspace(tree), _PDSamplingWorkspace()
    values = Vector{Float64}(undef, reps)
    for rep in 1:reps
        tips = _pd_sample!(sampler, probabilities, ntips, rng)
        values[rep] = _pd_site!(workspace, tree, tips, variant)
    end
    return values
end
