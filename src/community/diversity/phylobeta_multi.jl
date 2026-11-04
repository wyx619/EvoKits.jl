"""
    PhylobetaMultiResult

Multiple-site phylogenetic beta decomposition from [`phylobeta_multi`](@ref).
`total`, `turnover` and `nestedness` are scalars for the selected collection.
Owned `sample_labels`, `sample_indices`, `richness` and `pd` follow selection
order; indices refer to the original input rows. `pooled_pd` is the length of
the selected branch union, including paths to the original root.
"""
struct PhylobetaMultiResult
    sample_labels::Vector{String}
    sample_indices::Vector{Int}
    richness::Vector{Int}
    pd::Vector{Float64}
    pooled_pd::Float64
    total::Float64
    turnover::Float64
    nestedness::Float64
    index_family::Symbol
end

# Compensated sums keep small positive contributions; the error also bounds
# uncertainty in PD ordering and in differences between adjacent sorted PDs.
@inline function _phylobeta_multi_add(value, correction, increment)
    updated = value + increment
    error = value >= increment ? (value - updated) + increment :
        (increment - updated) + value
    return updated, correction + error, abs(error)
end

function _phylobeta_multi_components(a, m, x, difference, family)
    scale = max(a, m, x)
    scale > 0 || return (NaN, NaN, NaN)
    shared, low, high = a / scale, m / scale, x / scale
    denom = family === :sorensen ? 2shared + low + high : shared + low + high
    total = (low + high) / denom
    # Turnover has its own scale: a very large maximum unique length must not
    # erase both its shared and minimum-unique terms.
    turnover_scale = max(a, m)
    turnover_scale > 0 || return (Float64(total), NaN, NaN)
    shared, low = a / turnover_scale, m / turnover_scale
    numerator = family === :sorensen ? low : 2low
    turnover_denom = shared + numerator
    turnover = numerator / turnover_denom
    nestedness = (shared / turnover_denom) * ((difference / scale) / denom)
    return Float64(total), Float64(turnover), Float64(nestedness)
end

# For an edge present at q selected sites, A receives (q-1)*length.
# Sort sites by PD. For this edge let R=sum(rank-1) over its present sites.
# Its coefficients in the sums of maximum/minimum unique lengths are
# R-q(q-1)/2 and q(n-1)-R-q(q-1)/2. Both are nonnegative integers, so long
# shared stems cancel before multiplication by length, without subtracting PDs.
function _phylobeta_multi_calculate(branches, lengths, family, ::Type{T}) where {T}
    n = length(branches)
    values, errors = zeros(T, n), zeros(T, n)
    for (site, edges) in enumerate(branches)
        value = correction = uncertainty = zero(T)
        for edge in edges
            value, correction, error = _phylobeta_multi_add(value, correction, T(lengths[edge]))
            uncertainty += error
        end
        values[site] = value + correction
        errors[site] = 4uncertainty
    end
    all(isfinite, values) || return nothing
    order = sortperm(values)
    difference = difference_correction = zero(T)
    for k in 1:(n-1)
        left, right = order[k], order[k+1]
        gap = values[right] - values[left]
        if branches[left] == branches[right]
            gap = zero(T)
        elseif T === Float64 && errors[left] + errors[right] > 8eps(Float64) * gap
            return nothing
        end
        increment = T(Int128(k) * (n-k)) * gap
        difference, difference_correction, _ =
            _phylobeta_multi_add(difference, difference_correction, increment)
    end
    difference += difference_correction

    occupancy = zeros(Int, length(lengths))
    rank_sums = zeros(Int128, length(lengths))
    for (rank, site) in enumerate(order), edge in branches[site]
        occupancy[edge] += 1
        rank_sums[edge] += rank - 1
    end
    a = m = x = pooled = zero(T)
    ac = mc = xc = pc = zero(T)
    for edge in eachindex(lengths)
        q = occupancy[edge]
        q == 0 && continue
        paired = Int128(q) * (q-1) ÷ 2
        low = Int128(q) * (n-1) - rank_sums[edge] - paired
        high = rank_sums[edge] - paired
        weight = T(lengths[edge])
        a, ac, _ = _phylobeta_multi_add(a, ac, T(q-1) * weight)
        m, mc, _ = _phylobeta_multi_add(m, mc, T(low) * weight)
        x, xc, _ = _phylobeta_multi_add(x, xc, T(high) * weight)
        pooled, pc, _ = _phylobeta_multi_add(pooled, pc, weight)
    end
    a, m, x, pooled = a + ac, m + mc, x + xc, pooled + pc
    all(isfinite, (a, m, x, pooled, difference)) || return nothing
    if T === Float64
        largest = max(a, m, x, difference)
        # Use wider arithmetic when independent scaling would discard terms
        # near Float64 underflow. The final ratios still return Float64.
        any(v -> 0 < v < largest * 1e-280, (a, m, x, difference)) && return nothing
    end
    total, turnover, nestedness = _phylobeta_multi_components(a, m, x, difference, family)
    return (; pd=Float64.(values), pooled_pd=Float64(pooled), total, turnover, nestedness)
end

function _phylobeta_multi(branches, lengths, family)
    result = _phylobeta_multi_calculate(branches, lengths, family, Float64)
    result === nothing || return result
    # 2304 bits cover the entire binary64 exponent span plus Int128 rank
    # coefficients. This fallback resolves close PDs, overflow and subnormals;
    # it still uses occupancy/ranks and never enumerates site pairs.
    return setprecision(BigFloat, 2304) do
        _phylobeta_multi_calculate(branches, lengths, family, BigFloat)
    end
end

function _phylobeta_multi_result(branches, lengths, labels, indices, richness, family)
    result = _phylobeta_multi(branches, lengths, family)
    return PhylobetaMultiResult(labels[indices], copy(indices), copy(richness),
        result.pd, result.pooled_pd, result.total, result.turnover, result.nestedness, family)
end

"""
    phylobeta_multi(community, tree; index_family=:sorensen, samples=nothing,
                    tip_labels=nothing, sample_labels=nothing)
    phylobeta_multi(data::CommunityData, tree; index_family=:sorensen, samples=nothing)
    phylobeta_multi(context::CommunityContext; index_family=:sorensen, samples=nothing)

Decompose rooted phylogenetic beta diversity of a selected collection into
scalar total, turnover and nestedness components. Supports `:sorensen` and
`:jaccard`, dense/CSC matrices, `CommunityData` and reusable `CommunityContext`.
Finite positive abundance denotes presence. This is an overall multiple-site
index; with two sites it reduces to [`phylobeta`](@ref).

`samples` accepts a vector or tuple of distinct row indices or string labels;
omit it to use every row. At least two sites are required. Selection order and
original row indices are retained. Bool masks, repeated sites and mixed
index/label selectors are errors. All input values, labels and branch lengths
are validated, including unselected rows. Matrices need one column per tree
tip; `CommunityData` aligns labels and supplies zero columns for absent tips.

PD includes paths to the original tree root. `pooled_pd` counts only branches
present at selected sites; unobserved species and unselected sites do not add
branches. Context subsets recompute occupancy within the selection, while
reusing the stored branch sets. Neither input nor context is mutated.

With A = sum of selected PD minus pooled PD, and M/X the sums of minimum/maximum
unique branch length over unordered site pairs, Sørensen total is
(M+X)/(2A+M+X) and turnover M/(A+M); Jaccard total is (M+X)/(A+M+X) and turnover
2M/(A+2M). Nestedness uses the defining product, retaining small differences.

Empty and zero-PD sites are retained. All-zero PD gives three NaNs. One
positive-PD site among zero-PD sites gives total=1 and NaN turnover/nestedness.
Other zero denominators give NaN for the affected components. Zero-length
branches contribute no length. Tree branch lengths must be finite and
non-negative, with finite total length, as for other community metrics.

Branch occupancy and PD ranks replace pair enumeration. For m selected sites,
e tree edges and b selected site-branch memberships, the calculation after
branch collection costs O(b + e + m log m), with O(e+m) auxiliary storage.
No dense community, site-by-branch or site-pair matrix is built. Close PDs or
extreme scales use a wider arithmetic fallback without changing this layout.
"""
function phylobeta_multi(community::AbstractMatrix{<:Real}, tree::CompactTree;
    index_family::Symbol = :sorensen, samples = nothing,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    sample_labels = nothing,
)
    _validate_community_beta_family(index_family)
    _validate_community_branch_lengths(tree)
    labels = _validated_labels(sample_labels, size(community, 1), "sample")
    indices = _beta_multi_samples(samples, labels)
    length(indices) >= 2 || throw(ArgumentError("phylobeta_multi requires at least two selected sites"))
    rows = _community_tip_rows(community, tree, tip_labels)
    selected = view(rows, indices)
    branches = _community_branch_sets(selected, tree)
    return _phylobeta_multi_result(branches, tree.edge_length, labels, indices,
        length.(selected), index_family)
end

phylobeta_multi(community::AbstractMatrix, tree::CompactTree; kwargs...) =
    throw(ArgumentError("community matrix must contain real numbers without missing values"))
