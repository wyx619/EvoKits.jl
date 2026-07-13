function _community_distance_inputs(community::AbstractMatrix{<:Real}, tree::CompactTree, tip_labels)
    abundance = _community_abundance(community, tree, tip_labels)
    distances = pairwise_distance(tree)
    return abundance, distances
end

function _weighted_mean_or_nan(values, weights)
    keep = isfinite.(values) .& (weights .> 0)
    any(keep) || return NaN
    total = sum(weights[keep])
    total > 0 || return NaN
    return sum(values[keep] .* weights[keep]) / total
end

"""Mean pairwise phylogenetic distances between communities.

Rows are communities and columns are tree tips. The returned matrix follows
the input row order; diagonal entries are within-community mean distances.
Empty communities produce `NaN` entries.
"""
function comdist(community::AbstractMatrix{<:Real}, tree::CompactTree;
    abundance_weighted::Bool = false,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    abundance, distances = _community_distance_inputs(community, tree, tip_labels)
    n = size(abundance, 1)
    result = fill(NaN, n, n)
    for i in 1:n
        wi = abundance[i, :]
        for j in 1:i
            wj = abundance[j, :]
            if abundance_weighted
                weights = wi * wj'
            else
                weights = (wi .> 0) * (wj .> 0)'
            end
            total = sum(weights)
            total > 0 || continue
            result[i, j] = result[j, i] = sum(distances .* weights) / total
        end
    end
    return result
end

function _nearest_directional(abundance_from, abundance_to, distances;
    exclude_conspecifics::Bool,
    abundance_weighted::Bool,
)
    from = findall(>(0), abundance_from)
    to = findall(>(0), abundance_to)
    isempty(from) && return NaN
    isempty(to) && return NaN
    values = Float64[]
    weights = Float64[]
    for i in from
        candidates = copy(to)
        if exclude_conspecifics
            filter!(j -> j != i, candidates)
        end
        isempty(candidates) && continue
        nearest = minimum(distances[i, candidates])
        push!(values, nearest)
        push!(weights, abundance_weighted ? abundance_from[i] : 1.0)
    end
    isempty(values) && return NaN
    return _weighted_mean_or_nan(values, weights)
end

"""Mean nearest-taxon distances between each pair of communities."""
function comdistnt(community::AbstractMatrix{<:Real}, tree::CompactTree;
    abundance_weighted::Bool = false,
    exclude_conspecifics::Bool = false,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    abundance, distances = _community_distance_inputs(community, tree, tip_labels)
    n = size(abundance, 1)
    result = fill(NaN, n, n)
    for i in 1:n
        for j in 1:i
            if i == j
                value = _nearest_directional(abundance[i, :], abundance[j, :], distances;
                    exclude_conspecifics = exclude_conspecifics,
                    abundance_weighted = abundance_weighted)
            else
                forward = _nearest_directional(abundance[i, :], abundance[j, :], distances;
                    exclude_conspecifics = exclude_conspecifics,
                    abundance_weighted = abundance_weighted)
                reverse = _nearest_directional(abundance[j, :], abundance[i, :], distances;
                    exclude_conspecifics = exclude_conspecifics,
                    abundance_weighted = abundance_weighted)
                value = if isfinite(forward) && isfinite(reverse)
                    abundance_weighted ?
                        _weighted_mean_or_nan([forward, reverse], [sum(abundance[i, :]), sum(abundance[j, :])]) :
                        (forward + reverse) / 2
                elseif isfinite(forward)
                    forward
                else
                    reverse
                end
            end
            result[i, j] = result[j, i] = value
        end
    end
    return result
end

function _community_edge_presence(abundance::AbstractMatrix{<:Real}, tree::CompactTree)
    edge_presence = falses(size(abundance, 1), tree.nedges)
    paths = _tools_tip_paths(tree)
    for i in axes(abundance, 1)
        for tip in findall(>(0), abundance[i, :])
            edge_presence[i, paths[tip]] .= true
        end
    end
    return edge_presence
end

"""Unweighted UniFrac distances for rooted trees with branch lengths."""
function unifrac(community::AbstractMatrix{<:Real}, tree::CompactTree;
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    all(isfinite, tree.edge_length) && all(>=(0), tree.edge_length) ||
        throw(ArgumentError("unifrac requires finite, non-negative branch lengths"))
    abundance = _community_abundance(community, tree, tip_labels)
    active = _community_edge_presence(abundance, tree)
    n = size(abundance, 1)
    result = fill(NaN, n, n)
    for i in 1:n
        for j in 1:i
            union_length = 0.0
            unique_length = 0.0
            for edge in 1:tree.nedges
                ai, aj = active[i, edge], active[j, edge]
                (ai || aj) || continue
                len = tree.edge_length[edge]
                union_length += len
                ai != aj && (unique_length += len)
            end
            result[i, j] = result[j, i] = union_length > 0 ? unique_length / union_length : 0.0
        end
    end
    return result
end

"""PhyloSor similarity, `2 * shared PD / (PD₁ + PD₂)`."""
function phylosor(community::AbstractMatrix{<:Real}, tree::CompactTree;
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    abundance = _community_abundance(community, tree, tip_labels)
    active = _community_edge_presence(abundance, tree)
    edge_lengths = tree.edge_length
    n = size(abundance, 1)
    result = fill(NaN, n, n)
    pd_values = [sum(edge_lengths[active[i, :]]) for i in 1:n]
    for i in 1:n
        for j in 1:i
            denom = pd_values[i] + pd_values[j]
            result[i, j] = result[j, i] = denom > 0 ?
                2 * sum(edge_lengths[active[i, :] .& active[j, :]]) / denom : 0.0
        end
    end
    return result
end

Base.@kwdef struct RaoDResult
    Dkk::Vector{Float64} = Float64[]
    Dkl::Matrix{Float64} = zeros(Float64, 0, 0)
    H::Matrix{Float64} = zeros(Float64, 0, 0)
    total::Float64 = NaN
    alpha::Float64 = NaN
    beta::Float64 = NaN
    Fst::Float64 = NaN
end

"""Rao's quadratic entropy and alpha/beta partition for communities."""
function raoD(community::AbstractMatrix{<:Real}, tree::CompactTree;
    abundance_weighted::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    tip_depths = tree.dist_from_root[tree.tip_ids]
    maximum(abs.(tip_depths .- first(tip_depths))) <= 1e-8 ||
        throw(ArgumentError("raoD requires an ultrametric tree"))
    abundance, distances = _community_distance_inputs(community, tree, tip_labels)
    distances ./= 2.0
    n = size(abundance, 1)
    compositions = Matrix{Float64}(undef, n, tree.ntips)
    for i in 1:n
        weights = abundance_weighted ? copy(abundance[i, :]) : Float64.(abundance[i, :] .> 0)
        total = sum(weights)
        compositions[i, :] .= total > 0 ? weights ./ total : NaN
    end
    dkk = [sum(distances .* (compositions[i, :] * compositions[i, :]')) for i in 1:n]
    dkl = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:n
        dkl[i, j] = sum(distances .* (compositions[i, :] * compositions[j, :]'))
    end
    h = copy(dkl)
    for i in 1:n, j in 1:n
        h[i, j] = dkl[i, j] - (dkk[i] + dkk[j]) / 2
    end
    totals = vec(sum(abundance, dims = 2))
    total_abundance = sum(totals)
    pooled = total_abundance > 0 ? vec(sum(abundance, dims = 1)) ./ total_abundance : fill(NaN, tree.ntips)
    total = sum(distances .* (pooled * pooled'))
    sample_weights = total_abundance > 0 ? totals ./ total_abundance : fill(NaN, n)
    alpha = sum(dkk .* sample_weights)
    beta = total - alpha
    return RaoDResult(Dkk = dkk, Dkl = dkl, H = h, total = total, alpha = alpha, beta = beta, Fst = beta / total)
end
