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
    rows = _community_tip_rows(abundance, tree, nothing)
    _walk_community_branches(rows, tree) do site, edge
        edge_presence[site, edge] = true
    end
    return edge_presence
end

# Merge sorted edge sets and accumulate shared/unique lengths independently.
# Computing unique lengths by PD - shared can lose small terminal branches.
function _community_branch_components(left, right, lengths)
    a, b, c = 0.0, 0.0, 0.0
    i, j = 1, 1
    @inbounds while i <= length(left) && j <= length(right)
        le, re = left[i], right[j]
        if le == re
            a += lengths[le]
            i += 1
            j += 1
        elseif le < re
            b += lengths[le]
            i += 1
        else
            c += lengths[re]
            j += 1
        end
    end
    @inbounds for k in i:length(left)
        b += lengths[left[k]]
    end
    @inbounds for k in j:length(right)
        c += lengths[right[k]]
    end
    return a, b, c
end

function _community_pairwise(f, branches, lengths)
    n = length(branches)
    shared = n >= 128 ? _community_shared_lengths(branches, lengths) : nothing
    result = shared === nothing ? Matrix{Float64}(undef, n, n) : shared
    pd_values = shared === nothing ? Float64[] : diag(shared)
    for i in 1:n, j in 1:i
        a, b, c = shared === nothing ?
            _community_branch_components(branches[i], branches[j], lengths) :
            _community_components_from_shared(branches, lengths, shared[i, j], pd_values, i, j)
        result[i, j] = result[j, i] = f(a, b, c)
    end
    return result
end

# Full output already requires O(nsites^2) space. Sparse multiplication avoids
# merging every branch set for every pair; selected-pair calls keep the merge
# path and do not construct this matrix.
function _community_shared_lengths(branches, lengths)
    pointers = Vector{Int}(undef, length(branches) + 1)
    pointers[1] = 1
    for site in eachindex(branches)
        pointers[site + 1] = pointers[site] + length(branches[site])
    end
    indices = Vector{Int}(undef, last(pointers) - 1)
    for site in eachindex(branches)
        copyto!(indices, pointers[site], branches[site], 1, length(branches[site]))
    end
    incidence = SparseMatrixCSC(length(lengths), length(branches), pointers,
        indices, ones(Float64, length(indices)))
    weighted = copy(incidence)
    for k in eachindex(nonzeros(weighted))
        nonzeros(weighted)[k] = lengths[rowvals(weighted)[k]]
    end
    dropzeros!(weighted)
    return Matrix(transpose(incidence) * weighted)
end

@inline function _community_components_from_shared(branches, lengths, a, pd_values, i, j)
    b, c = pd_values[i] - a, pd_values[j] - a
    # Recompute unique branches directly near cancellation, including nested
    # or zero-PD sites and tiny terminal branches attached to a very long stem.
    if min(b, c) <= 1e-7 * max(pd_values[i], pd_values[j])
        return _community_branch_components(branches[i], branches[j], lengths)
    end
    return a, b, c
end

@inline function _community_scaled_components(a, b, c)
    scale = max(a, b, c)
    return scale > 0 ? (a / scale, b / scale, c / scale) : (0.0, 0.0, 0.0)
end

"""Unweighted UniFrac distances for rooted trees with branch lengths."""
function unifrac(community::AbstractMatrix{<:Real}, tree::CompactTree;
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    _validate_community_branch_lengths(tree)
    rows = _community_tip_rows(community, tree, tip_labels)
    branches = _community_branch_sets(rows, tree)
    return _unifrac(branches, tree.edge_length)
end

function _unifrac(branches, lengths)
    return _community_pairwise(branches, lengths) do a, b, c
        a, b, c = _community_scaled_components(a, b, c)
        denom = a + b + c
        denom > 0 ? (b + c) / denom : 0.0
    end
end

"""PhyloSor similarity, `2 * shared PD / (PD₁ + PD₂)`."""
function phylosor(community::AbstractMatrix{<:Real}, tree::CompactTree;
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    _validate_community_branch_lengths(tree)
    rows = _community_tip_rows(community, tree, tip_labels)
    branches = _community_branch_sets(rows, tree)
    return _phylosor(branches, tree.edge_length)
end

function _phylosor(branches, lengths)
    return _community_pairwise(branches, lengths) do a, b, c
        a, b, c = _community_scaled_components(a, b, c)
        denom = 2a + b + c
        denom > 0 ? 2a / denom : 0.0
    end
end

"""Pairwise beta decomposition; arrays are matrices, or vectors when pairs is supplied."""
struct PhylobetaResult{T<:AbstractArray{Float64}}
    sample_labels::Vector{String}
    total::T
    turnover::T
    nestedness::T
    index_family::Symbol
    pairs::Union{Nothing, Vector{Tuple{Int, Int}}}
end

@inline function _community_beta_components(a, b, c, family::Symbol)
    # Turnover only depends on shared length and the smaller unique length.
    # Scale these independently: a large other unique length can otherwise
    # underflow both to zero and turn a defined turnover into NaN.
    m = min(b, c)
    turnover_scale = max(a, m)
    shared, unique = turnover_scale > 0 ?
        (a / turnover_scale, m / turnover_scale) : (0.0, 0.0)
    turnover_denom = family === :sorensen ? shared + unique : shared + 2unique
    turnover = turnover_denom > 0 ?
        (family === :sorensen ? unique : 2unique) / turnover_denom : NaN
    a, b, c = _community_scaled_components(a, b, c)
    denom = family === :sorensen ? 2a + b + c : a + b + c
    total = denom > 0 ? (b + c) / denom : NaN
    nestedness = max(0.0, total - turnover)
    return total, turnover, nestedness
end

"""
    phylobeta(community, tree; index_family=:sorensen, pairs=nothing,
              tip_labels=nothing, sample_labels=nothing)

Partition rooted phylogenetic dissimilarity into turnover and nestedness.
Supports :sorensen and :jaccard. Returns three symmetric matrices in input
row order, or three vectors in the supplied order of (row_i, row_j) pairs.
The pairs option avoids constructing full site-by-site result matrices.

Zero denominators return NaN. Comparing an empty (or zero-PD) site with a
positive-PD site gives total=1 and undefined turnover/nestedness. Two zero-PD
sites give NaN for all components. Nonzero-PD diagonal entries are zero.
"""
function phylobeta(community::AbstractMatrix{<:Real}, tree::CompactTree;
    index_family::Symbol = :sorensen,
    pairs = nothing,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    _validate_community_beta_family(index_family)
    _validate_community_branch_lengths(tree)
    n = size(community, 1)
    labels = _tools_sample_labels(n, sample_labels)
    selected = _community_pairs(pairs, n)
    rows = _community_tip_rows(community, tree, tip_labels)
    branches = _community_branch_sets(rows, tree)
    total, turnover, nestedness = _phylobeta(branches, tree.edge_length, index_family, selected)
    return PhylobetaResult(labels, total, turnover, nestedness, index_family, selected)
end

function _validate_community_beta_family(family::Symbol)
    family in (:sorensen, :jaccard) ||
        throw(ArgumentError("index_family must be :sorensen or :jaccard"))
    return nothing
end

function _community_pairs(pairs, n)
    pairs === nothing && return nothing
    result = Tuple{Int, Int}[]
    for pair in pairs
        (pair isa Tuple || pair isa AbstractVector) && length(pair) == 2 &&
            all(x -> x isa Integer && !(x isa Bool) && 1 <= x <= n, pair) ||
            throw(ArgumentError("pairs must contain two valid integer row indices per pair"))
        push!(result, (Int(pair[1]), Int(pair[2])))
    end
    return result
end

function _phylobeta(branches, lengths, index_family, selected)
    n = length(branches)
    if selected === nothing
        shared = n >= 128 ? _community_shared_lengths(branches, lengths) : nothing
        total = shared === nothing ? Matrix{Float64}(undef, n, n) : shared
        turnover, nestedness = (Matrix{Float64}(undef, n, n) for _ in 1:2)
        pd_values = shared === nothing ? Float64[] : diag(shared)
        for i in 1:n, j in 1:i
            a, b, c = shared === nothing ?
                _community_branch_components(branches[i], branches[j], lengths) :
                _community_components_from_shared(branches, lengths, shared[i, j], pd_values, i, j)
            t, v, u = _community_beta_components(a, b, c, index_family)
            total[i, j] = total[j, i] = t
            turnover[i, j] = turnover[j, i] = v
            nestedness[i, j] = nestedness[j, i] = u
        end
    else
        total, turnover, nestedness = (Vector{Float64}(undef, length(selected)) for _ in 1:3)
        for (k, (i, j)) in enumerate(selected)
            a, b, c = _community_branch_components(branches[i], branches[j], lengths)
            total[k], turnover[k], nestedness[k] = _community_beta_components(a, b, c, index_family)
        end
    end
    return total, turnover, nestedness
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
