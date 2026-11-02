"""
    PhyloregionContext

Owned snapshot created by [`prepare_phyloregion`](@ref). Contains the original
packed distances, site labels, and complete hierarchy. Treat its arrays as
read-only; cuts and regionalization results own independent output arrays.
"""
struct PhyloregionContext
    sample_labels::Vector{String}
    distances::Vector{Float64}
    merge::Matrix{Int}
    heights::Vector{Float64}
    method::Symbol
    component::Symbol
end

function _prepare_phyloregion(distances, sample_labels, method, component)
    method = _region_method(method)
    n = size(distances, 1)
    labels = _community_data_labels(sample_labels, n, "sample")
    packed = _region_distance_buffer(distances)
    merge, heights = _region_hierarchy!(copy(packed), n, method)
    return PhyloregionContext(labels, packed, merge, heights, method, component)
end

"""
    prepare_phyloregion(distances; sample_labels=nothing, method=:average)
    prepare_phyloregion(beta; component=:total, method=:average)

Validate and copy complete distances, then build the hierarchy once. Reuse the
returned `PhyloregionContext` with `cutree(context; k)`, `phyloregion(context; k)`
and `optimal_phyloregion(context; ...)`. Inputs follow `phyloregion` rules.
The snapshot retains n*(n-1)/2 Float64 original distances plus O(n) tree/labels;
construction uses a second packed triangle as temporary clustering workspace.
Later cuts do not rebuild the tree. Regional summaries still visit site pairs.
"""
function prepare_phyloregion(distances::AbstractMatrix{<:Real};
    sample_labels = nothing, method::Symbol = :average,
)
    return _prepare_phyloregion(distances, sample_labels, method, :distance)
end

function prepare_phyloregion(beta::Union{BetaDissResult, PhylobetaResult};
    component::Symbol = :total, method::Symbol = :average,
)
    distances = _region_beta_distances(beta, component)
    return _prepare_phyloregion(distances, beta.sample_labels, method, component)
end

"""
    cutree(tree::Union{PhyloregionContext, PhyloregionResult}; k)

Cut an existing full hierarchy into `k` regions, returning an independent
membership vector in input-row order. Region IDs follow the first input site.
Cuts use merge count, including when centroid/median heights have inversions.
No distance matrix or regional summaries are allocated; work is O(n).
"""
function cutree(tree::Union{PhyloregionContext, PhyloregionResult}; k)
    n = length(tree.sample_labels)
    k = _region_validate_k(k, n)
    membership, _ = _region_cut(tree.merge, n, k)
    return membership
end

"""
    phyloregion(context::PhyloregionContext; k)

Reuse a prepared hierarchy and original distances for a new cut and regional
summaries. Does not recluster. Returns independently owned arrays, including
the full tree; pair summaries take O(n²) time and the region matrix O(k²) space.
"""
function phyloregion(context::PhyloregionContext; k)
    n = length(context.sample_labels)
    k = _region_validate_k(k, n)
    membership, sizes = _region_cut(context.merge, n, k)
    regions, distinctiveness = _region_summaries(context.distances, membership, k)
    return PhyloregionResult(copy(context.sample_labels), membership, sizes, regions,
        distinctiveness, copy(context.merge), copy(context.heights), k,
        context.method, context.component)
end
