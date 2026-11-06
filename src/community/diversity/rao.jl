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

# Pool on the original tree first, preserving original-abundance weights and
# tip order. Include that row in the existing connecting-tree mapping, then
# remove it from the site rows without modifying the prepared input records.
function _rao_prepare(rows, tree, pooled; strategy = :auto)
    strategy in (:auto, :full, :virtual) || throw(ArgumentError("unknown private RaoD strategy"))
    targets = count(!isempty, rows.nodes) + !isempty(pooled.nodes)
    virtual = !isempty(pooled.nodes) && (strategy === :virtual ||
        (strategy === :auto && _comdist_use_virtual(rows, tree, targets)))
    if virtual
        mapped, topology = _comdist_virtual((; nodes = [rows.nodes; [pooled.nodes]]), tree)
        pooled = merge(pooled, (; nodes = pop!(mapped.nodes)))
        return merge(rows, mapped), topology, pooled
    end
    isempty(rows.lengths) && (rows = merge(rows, (; lengths = _comdist_lengths(tree))))
    return rows, tree, pooled
end

# Evaluate only the diagonal of the existing distance kernel. Each site needs
# one target pass and its own records; no other source site is visited.
function _rao_diagonal(rows, tree, weights, workspace)
    dkk = fill(NaN, length(rows.nodes))
    for site in eachindex(rows.nodes)
        isempty(rows.nodes[site]) && continue
        _comdist_target!(workspace, tree, rows.nodes[site], weights[site], rows.lengths)
        value = 0.0
        for (node, weight) in zip(rows.nodes[site], weights[site])
            value += weight * (workspace.down[node] + workspace.outside_distance[node])
        end
        dkk[site] = value / 2
    end
    return dkk
end

"""
    raoD(community, tree; abundance_weighted=true, pairwise=true, tip_labels=nothing)
    raoD(data::CommunityData, tree; abundance_weighted=true, pairwise=true)

Rao's quadratic entropy and alpha/beta partition on an ultrametric tree,
using patristic distance / 2. Returns `RaoDResult` with per-site `Dkk`, site-pair
`Dkl`, `H[i,j] = Dkl[i,j] - (Dkk[i]+Dkk[j])/2`, pooled `total`, site-weighted
`alpha`, `beta = total-alpha`, and `Fst = beta/total`.
Set `pairwise=false` to compute only `Dkk` and the pooled partition, returning
independent 0×0 `Dkl` and `H` matrices without evaluating site pairs.

Site compositions use abundance by default, or equal weights over present
species with `abundance_weighted=false`. Both options pool the original
abundances and weight sites by their original abundance totals. Consequently
the incidence option can give negative beta/Fst. A zero pooled diversity gives
undefined Fst (NaN when beta is also zero).

Empty sites have NaN Dkk and propagate NaN to alpha/beta/Fst even with zero
site weight; total still uses any nonempty sites. In pairwise mode, their
Dkl/H rows and columns are NaN. All-empty input has NaN total; zero input
rows give alpha=0 and empty arrays. Labels follow the community input
conventions. Branch lengths must be finite and non-negative with a finite
sum; the ultrametric check retains the existing absolute root-to-tip
tolerance of 1e-8.

CSC input is traversed directly. Reuses the branch-distance kernel and one
tree workspace, without dense community, tip-distance or tip-pair weight
matrices. Sparse observed tip sets can use a connecting tree, contracting
unobserved/unary paths by summing original branches. Original-tree validation
and pooling precede compression; dense occupancy retains full-tree passes.
With the default `pairwise=true`, Dkl and H are full site-pair
matrices; `pairwise=false` uses linear auxiliary and output storage. Inputs
are not mutated.
"""
function raoD(community::AbstractMatrix{<:Real}, tree::CompactTree;
    abundance_weighted::Bool = true,
    pairwise::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    tip_depths = tree.dist_from_root[tree.tip_ids]
    maximum(abs.(tip_depths .- first(tip_depths))) <= 1e-8 ||
        throw(ArgumentError("raoD requires an ultrametric tree"))
    rows = _community_distance_rows(community, tree, tip_labels, true; need_lengths = false)
    return _rao_result(rows, tree; abundance_weighted, pairwise)
end

# The private strategy override supports comparisons of both topologies; it
# does not add a public keyword or alter either result mode's contract.
function _rao_result(rows, tree; abundance_weighted = true, pairwise = true, strategy = :auto)
    pooled = _rao_pool(rows, tree)
    rows, tree, pooled = _rao_prepare(rows, tree, pooled; strategy)
    weights = abundance_weighted ? rows.weights :
        [isempty(nodes) ? Float64[] : fill(1.0 / length(nodes), length(nodes)) for nodes in rows.nodes]
    workspace = _ComdistWorkspace(tree.nnodes)
    if pairwise
        dkl = _comdist(rows, tree, weights, workspace)
        dkl ./= 2
        dkk = diag(dkl)
        n = length(dkk)
        h = Matrix{Float64}(undef, n, n)
        for i in 1:n, j in 1:i
            h[i, j] = h[j, i] = dkl[i, j] - (dkk[i] + dkk[j]) / 2
        end
    else
        dkk = _rao_diagonal(rows, tree, weights, workspace)
        dkl, h = zeros(0, 0), zeros(0, 0)
    end
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
