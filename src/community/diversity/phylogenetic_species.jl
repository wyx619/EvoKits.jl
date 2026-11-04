Base.@kwdef struct PSVResult
    psv::Vector{Float64} = Float64[]
    richness::Vector{Int} = Int[]
    variance::Vector{Float64} = Float64[]
    scale_vcv::Bool = true
end

# PSV is the mean over unordered pairs of (Cii+Cjj)/2-Cij: distance/2
# for raw Brownian covariance, or 1-Cij for correlation covariance.
# Deficits 1-sqrt(h_parent/h_tip) preserve correlations close to one.
function _psv_geometry(tree, scaled)
    lengths = zeros(tree.nnodes)
    for edge in 1:tree.nedges
        lengths[tree.child_of_edge[edge]] = tree.edge_length[edge]
    end
    factors, offsets = ones(tree.nnodes), copy(lengths)
    if scaled
        heights = zeros(tree.nnodes)
        for node in tree.preorder
            parent = tree.parent_of_node[node]
            parent == 0 && continue
            heights[node] = heights[parent] + lengths[node]
        end
        all(node -> heights[node] > 0, tree.tip_ids) ||
            throw(ArgumentError("tree has non-positive root-to-tip variance"))
        for node in tree.preorder
            parent = tree.parent_of_node[node]
            parent == 0 && continue
            if heights[node] > 0
                factors[node] = sqrt(heights[parent] / heights[node])
                offsets[node] = (lengths[node] / heights[node]) / (1 + factors[node])
            else
                factors[node], offsets[node] = 1.0, 0.0
            end
        end
    end
    return (; factors, offsets, scaled)
end

# Merge populations using SD rather than an unnormalized sum of squares,
# avoiding overflow from squared distances or the number of tip pairs.
@inline function _psv_merge(na, ma, sa, nb, mb, sb)
    na == 0 && return (nb, mb, sb)
    nb == 0 && return (na, ma, sa)
    n = na + nb
    a, b = na / n, nb / n
    delta = mb - ma
    return n, ma + b * delta, hypot(sqrt(a) * sa, sqrt(b) * sb, sqrt(a * b) * delta)
end

@inline function _psv_merge_mean(na, ma, nb, mb)
    na == 0 && return (nb, mb)
    nb == 0 && return (na, ma)
    return na + nb, ma + (nb / (na + nb)) * (mb - ma)
end

function _psv_moments!(tree, tips, geometry, counts, means, sds; central = false)
    fill!(counts, 0)
    fill!(means, 0.0)
    central && fill!(sds, 0.0)
    for tip in tips
        counts[tree.tip_ids[tip]] = 1
    end
    pairs, pair_mean, pair_sd = 0, 0.0, 0.0
    for node in tree.postorder
        n, m, s = counts[node], 0.0, 0.0
        for child in tree.children[node]
            nc = counts[child]
            nc == 0 && continue
            mc = geometry.offsets[child] + geometry.factors[child] * means[child]
            cross_mean = geometry.scaled ? m + (1 - m) * mc : m + mc
            if central
                sc = geometry.factors[child] * sds[child]
                cross_sd = geometry.scaled ? hypot((1 - mc) * s, (1 - m) * sc, s * sc) : hypot(s, sc)
                pairs, pair_mean, pair_sd = _psv_merge(pairs, pair_mean, pair_sd,
                    n * nc, cross_mean, cross_sd)
                n, m, s = _psv_merge(n, m, s, nc, mc, sc)
            else
                pairs, pair_mean = _psv_merge_mean(pairs, pair_mean, n * nc, cross_mean)
                n, m = _psv_merge_mean(n, m, nc, mc)
            end
        end
        counts[node], means[node] = n, m
        central && (sds[node] = s)
    end
    return pair_mean, pair_sd
end

# Mean contrast from each tip to all other tips. Prefix/suffix sibling means
# avoid subtracting two nearly equal subtree sums on an unbalanced tree.
function _psv_row_sd(tree, geometry, counts, means, pair_mean)
    outside = zeros(tree.nnodes)
    suffix_count, suffix_mean = zeros(Int, tree.nnodes), zeros(tree.nnodes)
    p = tree.ntips
    for parent in tree.preorder
        after_count, after_mean = 0, 0.0
        for child in Iterators.reverse(tree.children[parent])
            suffix_count[child], suffix_mean[child] = after_count, after_mean
            child_mean = geometry.offsets[child] + geometry.factors[child] * means[child]
            after_count, after_mean = _psv_merge_mean(after_count, after_mean, counts[child], child_mean)
        end
        before_count, before_mean = 0, 0.0
        for child in tree.children[parent]
            n, m = _psv_merge_mean(p - counts[parent], outside[parent], before_count, before_mean)
            n, m = _psv_merge_mean(n, m, suffix_count[child], suffix_mean[child])
            n > 0 && (outside[child] = geometry.offsets[child] + geometry.factors[child] * m)
            child_mean = geometry.offsets[child] + geometry.factors[child] * means[child]
            before_count, before_mean = _psv_merge_mean(before_count, before_mean, counts[child], child_mean)
        end
    end
    norm = 0.0
    normalization = sqrt(p)
    for tip in tree.tip_ids
        norm = hypot(norm, (outside[tip] - pair_mean) / normalization)
    end
    return norm
end

# Exact finite-population variance of the mean pair contrast at richness m.
# Nonnegative centered pair/row coefficients replace picante's asymmetric SS2.
function _psv_variances(tree, geometry, counts, means, sds)
    p = tree.ntips
    result = fill(NaN, p)
    p <= 1 && return result
    pair_mean, pair_sd = _psv_moments!(tree, 1:p, geometry, counts, means, sds; central = true)
    row_sd = _psv_row_sd(tree, geometry, counts, means, pair_mean)
    conversion = geometry.scaled ? 1.0 : 0.5
    result[p] = 0.0
    if p == 3
        result[2] = (conversion * pair_sd)^2
    elseif p >= 4
        for m in 2:(p - 1)
            pairs = m * (m - 1) / 2
            denominator = pairs * (p - 2) * (p - 3)
            a = (p - m) * (p - m - 1) / denominator
            b = 2.0 * (p - 1) * (m - 2) * (p - m) / denominator
            result[m] = (conversion * sqrt(a) * pair_sd)^2 + (conversion * sqrt(b) * row_sd)^2
        end
    end
    return result
end

"""
    psv(community, tree; scale_vcv=true, compute_var=true, tip_labels=nothing)

Phylogenetic species variability in input site order. Positive values denote
presence. Uses Brownian correlation by default, or raw covariance with
`scale_vcv=false`, without constructing VCV or site submatrices. Empty and
singleton sites give NaN. Scaled calls require positive root-to-tip lengths.

Variance is exact for uniform sampling without replacement at each site's
richness from ALL tree tips, including unobserved tips. A full pool has zero
variance. This fixes the previous leaf-order-dependent formula; raw covariance
variance also includes unequal diagonal contributions. `compute_var=false`
leaves variance as NaN. Inputs are not mutated.
"""
function psv(community::AbstractMatrix{<:Real}, tree::CompactTree;
    scale_vcv::Bool = true,
    compute_var::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    _validate_community_branch_lengths(tree)
    rows = _community_tip_rows(community, tree, tip_labels)
    geometry = _psv_geometry(tree, scale_vcv)
    counts, means, sds = zeros(Int, tree.nnodes), zeros(tree.nnodes), zeros(tree.nnodes)
    variance_by_richness = compute_var ? _psv_variances(tree, geometry, counts, means, sds) : Float64[]
    values, variances = fill(NaN, length(rows)), fill(NaN, length(rows))
    richness = length.(rows)
    for (site, tips) in enumerate(rows)
        length(tips) > 1 || continue
        value, _ = _psv_moments!(tree, tips, geometry, counts, means, sds)
        values[site] = scale_vcv ? value : value / 2
        compute_var && (variances[site] = variance_by_richness[length(tips)])
    end
    return PSVResult(psv = values, richness = richness, variance = variances, scale_vcv = scale_vcv)
end

function psr(community::AbstractMatrix{<:Real}, tree::CompactTree; kwargs...)
    result = psv(community, tree; kwargs...)
    return (psr = result.psv .* result.richness, richness = result.richness,
        variance = result.variance .* result.richness .^ 2)
end
