Base.@kwdef struct RaoDResult
    Dkk::Vector{Float64} = Float64[]
    Dkl::Matrix{Float64} = zeros(Float64, 0, 0)
    H::Matrix{Float64} = zeros(Float64, 0, 0)
    total::Float64 = NaN
    alpha::Float64 = NaN
    beta::Float64 = NaN
    Fst::Float64 = NaN
end

# Pool original abundances even when site compositions use incidence. Scale
# before summing so finite counts need not have a representable raw total.
function _rao_pool(rows, tree)
    site_weights = zeros(length(rows.nodes))
    pooled = zeros(tree.nnodes)
    scale = maximum(rows.scales; init = 0)
    scale > 0 || return (; nodes = Int32[], weights = Float64[], site_weights = fill(NaN, length(site_weights)))
    for site in eachindex(rows.nodes)
        for (node, value) in zip(rows.nodes[site], rows.raw[site])
            weight = _community_distance_ratio(value, scale)
            pooled[node] += weight
            site_weights[site] += weight
        end
    end
    total = sum(site_weights)
    site_weights ./= total
    nodes = Int32[node for node in tree.tip_ids if pooled[node] > 0]
    weights = Float64[pooled[node] / total for node in nodes]
    return (; nodes, weights, site_weights)
end

"""
    raoD(community, tree; abundance_weighted=true, tip_labels=nothing)
    raoD(data::CommunityData, tree; abundance_weighted=true)

Rao's quadratic entropy and alpha/beta partition on an ultrametric tree,
using patristic distance / 2. Returns `RaoDResult` with per-site `Dkk`, site-pair
`Dkl`, `H[i,j] = Dkl[i,j] - (Dkk[i]+Dkk[j])/2`, pooled `total`, site-weighted
`alpha`, `beta = total-alpha`, and `Fst = beta/total`.

Site compositions use abundance by default, or equal weights over present
species with `abundance_weighted=false`. Both options pool the original
abundances and weight sites by their original abundance totals. Consequently
the incidence option can give negative beta/Fst. A zero pooled diversity gives
undefined Fst (NaN when beta is also zero).

Empty sites have NaN Dkk and NaN Dkl/H rows and columns. They also propagate
NaN to alpha/beta/Fst even with zero site weight; total still uses any nonempty
sites. All-empty input has NaN total; zero input rows give alpha=0 and empty
arrays. Labels follow the community input conventions. Branch lengths must
be finite and non-negative with a finite sum; the ultrametric check retains
the existing absolute root-to-tip tolerance of 1e-8.

CSC input is traversed directly. Reuses the branch-distance kernel and one
tree workspace, without dense community, tip-distance or tip-pair weight
matrices. Dkl and H remain full site-pair matrices. Inputs are not mutated.
"""
function raoD(community::AbstractMatrix{<:Real}, tree::CompactTree;
    abundance_weighted::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    tip_depths = tree.dist_from_root[tree.tip_ids]
    maximum(abs.(tip_depths .- first(tip_depths))) <= 1e-8 ||
        throw(ArgumentError("raoD requires an ultrametric tree"))
    rows = _community_distance_rows(community, tree, tip_labels, true)
    weights = abundance_weighted ? rows.weights :
        [isempty(nodes) ? Float64[] : fill(1.0 / length(nodes), length(nodes)) for nodes in rows.nodes]
    workspace = _ComdistWorkspace(tree.nnodes)
    dkl = _comdist(rows, tree, weights, workspace)
    dkl ./= 2
    dkk = diag(dkl)
    n = length(dkk)
    h = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:i
        h[i, j] = h[j, i] = dkl[i, j] - (dkk[i] + dkk[j]) / 2
    end
    pooled = _rao_pool(rows, tree)
    total = NaN
    if !isempty(pooled.nodes)
        _comdist_target!(workspace, tree, pooled.nodes, pooled.weights, rows.lengths)
        total = 0.0
        for (node, weight) in zip(pooled.nodes, pooled.weights)
            total += weight * (workspace.down[node] + workspace.outside_distance[node])
        end
        total /= 2
    end
    alpha = sum(dkk .* pooled.site_weights)
    beta = total - alpha
    return RaoDResult(Dkk = dkk, Dkl = dkl, H = h, total = total,
        alpha = alpha, beta = beta, Fst = beta / total)
end
