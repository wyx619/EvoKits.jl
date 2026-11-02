"""
    PhyloregionResult

Regionalization returned by [`phyloregion`](@ref). `membership` assigns
each input site to regions `1:k` in first-site order; `cluster_sizes` follows
that region order. `region_distances` has a zero diagonal and mean cross-site
distances off diagonal. `distinctiveness` is each region's unweighted mean
distance to the other regions (NaN for k=1).

`merge` and `heights` describe the full hierarchy, including merges above the
requested cut: negative merge entries identify input rows, and positive entries
refer to earlier merge rows, as in R hclust. Arrays are owned by the result.
`component` is :distance for matrix input, or the selected beta component.
"""
struct PhyloregionResult
    sample_labels::Vector{String}
    membership::Vector{Int}
    cluster_sizes::Vector{Int}
    region_distances::Matrix{Float64}
    distinctiveness::Vector{Float64}
    merge::Matrix{Int}
    heights::Vector{Float64}
    k::Int
    method::Symbol
    component::Symbol
end

const _REGION_METHODS = (:average, :single, :complete, :ward_D, :ward_D2,
    :mcquitty, :median, :centroid)

function _region_method(method::Symbol)
    canonical = method === Symbol("ward.D") ? :ward_D :
        method === Symbol("ward.D2") ? :ward_D2 : method
    canonical in _REGION_METHODS || throw(ArgumentError("unsupported linkage method: $method"))
    return canonical
end

# Packed upper triangle, one Float64 per unordered pair. Physical cluster slots
# retain their smallest input-row index throughout nearest-neighbor chaining.
@inline function _region_pair_index(a::Int, b::Int)
    lo, hi = minmax(a, b)
    return (hi - 1) * (hi - 2) ÷ 2 + lo
end

function _region_distance_buffer(distances)
    Base.require_one_based_indexing(distances)
    n = size(distances, 1)
    size(distances, 2) == n || throw(ArgumentError("distance matrix must be square"))
    n > 0 || throw(ArgumentError("regionalization requires at least one site"))
    packed = Vector{Float64}(undef, n * (n - 1) ÷ 2)
    for j in 1:n
        distances[j, j] == 0 || throw(ArgumentError("distance diagonal must be zero and finite; remove undefined sites explicitly"))
        for i in 1:(j - 1)
            value, reverse = distances[i, j], distances[j, i]
            isfinite(value) && isfinite(reverse) && value >= 0 && reverse >= 0 ||
                throw(ArgumentError("distances must be finite and non-negative; remove undefined sites explicitly"))
            value == reverse || throw(ArgumentError("distance matrix must be symmetric"))
            converted = Float64(value)
            isfinite(converted) && (value == 0 || converted > 0) ||
                throw(ArgumentError("distances must be representable as finite Float64 values without underflow to zero"))
            packed[_region_pair_index(i, j)] = converted
        end
    end
    return packed
end

# Convex interpolation of nonnegative finite distances without a weighted sum
# that can overflow. Retain the bounds even at the Float64 rounding limits.
@inline function _region_average(left, right, left_size, right_size)
    lo, hi = minmax(left, right)
    hi_weight = left <= right ? right_size : left_size
    return clamp(lo + (hi - lo) * (hi_weight / (left_size + right_size)), lo, hi)
end

@inline _region_update(x, y, d, a, b, c, ::Val{:average}) = _region_average(x, y, a, b)
@inline _region_update(x, y, d, a, b, c, ::Val{:single}) = min(x, y)
@inline _region_update(x, y, d, a, b, c, ::Val{:complete}) = max(x, y)
@inline _region_update(x, y, d, a, b, c, ::Val{:mcquitty}) = _region_average(x, y, 1, 1)

# Scale each update independently. In particular, Ward D2 never squares the
# unscaled input, which could overflow or underflow for valid finite distances.
@inline function _region_update(x, y, d, a, b, c, ::Val{M}) where {M}
    scale = max(x, y, d)
    scale == 0 && return 0.0
    u, v, w = x / scale, y / scale, d / scale
    if M === :ward_D || M === :ward_D2
        p, q, r = (a + c) / (a + b + c), (b + c) / (a + b + c), c / (a + b + c)
        value = M === :ward_D2 ? p * u^2 + q * v^2 - r * w^2 : p * u + q * v - r * w
    elseif M === :centroid
        p, q = a / (a + b), b / (a + b)
        value = p * u + q * v - p * q * w
    else # :median, R's WPGMC recurrence on the supplied dissimilarities.
        value = (u + v) / 2 - w / 4
    end
    value >= -64eps(Float64) || throw(ArgumentError("linkage update produced a negative distance"))
    value = max(0.0, value)
    updated = (M === :ward_D2 ? sqrt(value) : value) * scale
    isfinite(updated) || throw(ArgumentError("linkage distances exceed finite Float64 range; rescale the input"))
    return M in (:ward_D, :ward_D2) ? max(updated, min(x, y)) : updated
end

function _region_chain!(packed, n, criterion::Val)
    active = trues(n)
    sizes = ones(Int, n)
    nodes = -collect(1:n)
    chain = Int[]
    raw_merge = Matrix{Int}(undef, n - 1, 2)
    raw_heights = Vector{Float64}(undef, n - 1)
    for step in 1:(n - 1)
        if isempty(chain)
            push!(chain, findfirst(active))
        end
        while true
            current = last(chain)
            previous = length(chain) > 1 ? chain[end - 1] : 0
            nearest = previous
            best = previous == 0 ? Inf : packed[_region_pair_index(current, previous)]
            for candidate in 1:n
                (!active[candidate] || candidate == current) && continue
                value = packed[_region_pair_index(current, candidate)]
                # A tie retains the predecessor, preventing equal-distance
                # cycles. With no predecessor, scan order chooses the lowest
                # active slot. Strictly closer candidates replace either.
                if value < best
                    nearest, best = candidate, value
                end
            end
            if nearest == previous
                pop!(chain)
                pop!(chain)
                keep, drop = minmax(current, nearest)
                raw_merge[step, 1], raw_merge[step, 2] = nodes[keep], nodes[drop]
                raw_heights[step] = best
                for other in 1:n
                    (!active[other] || other == keep || other == drop) && continue
                    first = _region_pair_index(keep, other)
                    second = _region_pair_index(drop, other)
                    packed[first] = _region_update(packed[first], packed[second], best,
                        sizes[keep], sizes[drop], sizes[other], criterion)
                end
                sizes[keep] += sizes[drop]
                active[drop] = false
                nodes[keep] = step
                break
            end
            push!(chain, nearest)
        end
    end
    # Reciprocal nearest neighbors can be found out of global height order.
    # These six methods are reducible: a stable height sort gives the hierarchy;
    # equal-height children precede parents in their raw creation order.
    order = sortperm(raw_heights; alg = Base.Sort.MergeSort)
    indices = invperm(order)
    merge = similar(raw_merge)
    for (step, raw) in enumerate(order), child in 1:2
        node = raw_merge[raw, child]
        merge[step, child] = node < 0 ? node : indices[node]
    end
    return merge, raw_heights[order]
end

# Centroid/median may have inversions: use global minimum selection, retain
# chronological merge order, and never sort their heights. Cache row minima
# to avoid scanning the entire triangle at every merge. Worst case O(n^3).
function _region_global!(packed, n, criterion::Val)
    active, sizes, nodes = trues(n), ones(Int, n), -collect(1:n)
    nearest, minima = zeros(Int, n), fill(Inf, n)
    merge, heights = Matrix{Int}(undef, n - 1, 2), Vector{Float64}(undef, n - 1)
    function refresh!(i)
        nearest[i], minima[i] = 0, Inf
        for j in (i + 1):n
            active[j] || continue
            value = packed[_region_pair_index(i, j)]
            if value < minima[i]
                nearest[i], minima[i] = j, value
            end
        end
    end
    for i in 1:n
        refresh!(i)
    end
    for step in 1:(n - 1)
        keep, best = 0, Inf
        for i in 1:n
            if active[i] && minima[i] < best
                keep, best = i, minima[i]
            end
        end
        drop = nearest[keep]
        merge[step, 1], merge[step, 2] = nodes[keep], nodes[drop]
        heights[step] = best
        active[drop] = false
        for other in 1:n
            (!active[other] || other == keep) && continue
            pair = _region_pair_index(keep, other)
            packed[pair] = _region_update(packed[pair], packed[_region_pair_index(drop, other)],
                best, sizes[keep], sizes[drop], sizes[other], criterion)
        end
        sizes[keep] += sizes[drop]
        nodes[keep] = step
        for i in 1:n
            active[i] || continue
            if i == keep || nearest[i] == keep || nearest[i] == drop
                refresh!(i)
            elseif i < keep
                value = packed[_region_pair_index(i, keep)]
                if value < minima[i] || (value == minima[i] && keep < nearest[i])
                    nearest[i], minima[i] = keep, value
                end
            end
        end
    end
    return merge, heights
end

_region_hierarchy!(packed, n, method) = method in (:centroid, :median) ?
    _region_global!(packed, n, Val(method)) : _region_chain!(packed, n, Val(method))

function _region_find_root!(parents, node)
    root = node
    while parents[root] != root
        root = parents[root]
    end
    while parents[node] != node
        following = parents[node]
        parents[node] = root
        node = following
    end
    return root
end

function _region_cut(merge, n, k)
    parents = collect(1:n)
    representatives = Vector{Int}(undef, n - 1)
    for step in 1:(n - k)
        a, b = merge[step, 1], merge[step, 2]
        left = _region_find_root!(parents, a < 0 ? -a : representatives[a])
        right = _region_find_root!(parents, b < 0 ? -b : representatives[b])
        keep, drop = minmax(left, right)
        parents[drop] = keep
        representatives[step] = keep
    end
    membership = Vector{Int}(undef, n)
    region_of_root = zeros(Int, n)
    sizes = Int[]
    for site in 1:n
        root = _region_find_root!(parents, site)
        if region_of_root[root] == 0
            push!(sizes, 0)
            region_of_root[root] = length(sizes)
        end
        region = region_of_root[root]
        membership[site] = region
        sizes[region] += 1
    end
    return membership, sizes
end

@inline _region_value(distances::AbstractMatrix, i, j) = Float64(distances[i, j])
@inline _region_value(distances::AbstractVector, i, j) = distances[_region_pair_index(i, j)]

function _region_summaries(distances, membership, k)
    regions = zeros(Float64, k, k)
    counts = zeros(Int, k * (k - 1) ÷ 2)
    for j in eachindex(membership), i in 1:(j - 1)
        a, b = minmax(membership[i], membership[j])
        a == b && continue
        pair = _region_pair_index(a, b)
        counts[pair] += 1
        # Online means avoid overflow when many large finite distances are
        # averaged, both within a region pair and across other regions.
        regions[a, b] += (_region_value(distances, i, j) - regions[a, b]) / counts[pair]
    end
    for b in 1:k, a in 1:(b - 1)
        regions[b, a] = regions[a, b]
    end
    distinctiveness = fill(NaN, k)
    if k > 1
        for a in 1:k
            value, count = 0.0, 0
            for b in 1:k
                a == b && continue
                count += 1
                value += (regions[a, b] - value) / count
            end
            distinctiveness[a] = value
        end
    end
    return regions, distinctiveness
end

function _region_validate_k(k, n)
    k isa Integer && !(k isa Bool) || throw(ArgumentError("k must be an integer"))
    1 <= k <= n || throw(ArgumentError("k must be between 1 and the number of sites"))
    return Int(k)
end

function _phyloregion(distances, k, sample_labels, method, component)
    method = _region_method(method)
    n = size(distances, 1)
    k = _region_validate_k(k, n)
    labels = _community_data_labels(sample_labels, n, "sample")
    packed = _region_distance_buffer(distances)
    merge, heights = _region_hierarchy!(packed, n, method)
    membership, sizes = _region_cut(merge, n, Int(k))
    regions, distinctiveness = _region_summaries(distances, membership, Int(k))
    return PhyloregionResult(labels, membership, sizes, regions, distinctiveness,
        merge, heights, Int(k), method, component)
end

"""
    phyloregion(distances; k, sample_labels=nothing, method=:average)
    phyloregion(beta::Union{BetaDissResult, PhylobetaResult}; k,
                component=:total, method=:average)

Partition sites by hierarchical clustering at an explicit number of regions
`k`, and summarize regional distances. Methods are `:average` (UPGMA, default),
`:single`, `:complete`, `:mcquitty` (WPGMA), `:ward_D`, `:ward_D2`, `:centroid`
and `:median`. `Symbol("ward.D")` / `Symbol("ward.D2")` are aliases normalized
to `:ward_D` / `:ward_D2`. Ward D2 squares distances internally and reports
square-root heights; other methods apply R hclust recurrences directly to the
supplied values. No automatic squaring is done for centroid/median. Ward's
minimum-variance interpretation requires Euclidean distances.
Matrix input must be nonempty,
square, exactly symmetric, non-negative, finite and have a zero diagonal.
Distances use Float64 and must not overflow or underflow to zero on conversion.
No Euclidean or triangle-inequality assumption is required. Labels must be
unique nonempty strings or integer identifiers, and retain input site order.

Beta input must contain full matrices (`pairs=nothing`); select `:total`,
`:turnover`, or `:nestedness`. Its site labels are retained. Undefined distances,
including NaN diagonals from empty sites, are errors: explicitly select a
defined subset before calling; no site is dropped and no value is imputed.

Returns `PhyloregionResult` with site `membership`, `cluster_sizes`, mean
cross-site `region_distances`, and each region's mean distance to the other
regions (`distinctiveness`). Regions receive equal weight in distinctiveness,
regardless of their site counts. Region diagonals are zero; k=1 gives NaN
distinctiveness. The result also contains the complete `merge`/`heights` tree.

The six monotone methods use an O(n^2) nearest-neighbor-chain algorithm with
one packed Float64 work buffer of n*(n-1)/2 distances. A fresh chain starts at
the active cluster containing
the earliest input site; equal nearest distances retain the predecessor if
present, otherwise the earliest active cluster wins. Equal-height merges keep
discovery order. This deterministic tie rule can differ from R hclust; ties
may also change under input-row permutation. Centroid/median use global
minimum selection with cached row minima (worst-case O(n^3)), choosing the
lexicographically first pair of active slots on a tie. Their heights may
decrease: chronological merge order is preserved and never height-sorted.
No tie tolerance is applied. A linkage update exceeding finite Float64 range
is an error; rescale input distances in that case.
Matrix/label inputs are not mutated and results own their arrays. Full input
and region output still require quadratic memory; selected pairs alone are
insufficient. Use `prepare_phyloregion` to reuse a tree, `select_linkage` to
compare methods, and `optimal_phyloregion` to evaluate cuts. Plotting is separate.
"""
function phyloregion(distances::AbstractMatrix{<:Real}; k,
    sample_labels = nothing, method::Symbol = :average,
)
    return _phyloregion(distances, k, sample_labels, method, :distance)
end

function phyloregion(beta::Union{BetaDissResult, PhylobetaResult}; k,
    component::Symbol = :total, method::Symbol = :average,
)
    distances = _region_beta_distances(beta, component)
    return _phyloregion(distances, k, beta.sample_labels, method, component)
end

function _region_beta_distances(beta, component)
    component in (:total, :turnover, :nestedness) ||
        throw(ArgumentError("component must be :total, :turnover or :nestedness"))
    beta.pairs === nothing || throw(ArgumentError("regionalization requires full beta matrices, not selected pairs"))
    distances = getproperty(beta, component)
    distances isa AbstractMatrix || throw(ArgumentError("regionalization requires a distance matrix"))
    return distances
end
