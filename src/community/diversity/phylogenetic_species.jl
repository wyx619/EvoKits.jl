Base.@kwdef struct PSVResult
    psv::Vector{Float64} = Float64[]
    richness::Vector{Int} = Int[]
    variance::Vector{Float64} = Float64[]
    scale_vcv::Bool = true
end

# PSV is the mean over unordered pairs of (Cii+Cjj)/2-Cij: distance/2
# for raw Brownian covariance, or 1-Cij for correlation covariance.
# Deficits 1-sqrt(h_parent/h_tip) preserve correlations close to one.
function _psv_geometry(tree, scaled; need_heights = false)
    lengths = zeros(tree.nnodes)
    for edge in 1:tree.nedges
        lengths[tree.child_of_edge[edge]] = tree.edge_length[edge]
    end
    factors, offsets = ones(tree.nnodes), copy(lengths)
    heights = scaled || need_heights ? zeros(tree.nnodes) : Float64[]
    if !isempty(heights)
        for node in tree.preorder
            parent = tree.parent_of_node[node]
            parent == 0 && continue
            heights[node] = heights[parent] + lengths[node]
        end
    end
    if scaled
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
    return (; factors, offsets, down_offsets = offsets, scaled, heights)
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
Sparse sites reuse observed connecting trees and site workspaces; dense sites
keep sequential tree traversal. PSR uses the same preparation through PSV.

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
    plan = _species_plan(rows, tree, geometry)
    if plan.index === nothing && plan.tree isa CompactTree
        # Keep the existing dense route, sharing its full-pool variance scratch
        # with site moments instead of allocating a second full workspace.
        counts, means, sds = zeros(Int, tree.nnodes), zeros(tree.nnodes), zeros(tree.nnodes)
        return _psv_values(rows, tree, geometry, counts, means, sds; compute_var)
    end
    # Compression serves site means only. Unobserved and singleton-only tips
    # remain in the original finite population used to calculate variance.
    variance_by_richness = compute_var ? _psv_variances(tree, geometry,
        zeros(Int, tree.nnodes), zeros(tree.nnodes), zeros(tree.nnodes)) : Float64[]
    return _species_psv_values(rows, plan, variance_by_richness)
end

function _psv_values(rows, tree, geometry, counts, means, sds; compute_var)
    variance_by_richness = compute_var ? _psv_variances(tree, geometry, counts, means, sds) : Float64[]
    values, variances = fill(NaN, length(rows)), fill(NaN, length(rows))
    richness = length.(rows)
    for (site, tips) in enumerate(rows)
        length(tips) > 1 || continue
        value, _ = _psv_moments!(tree, tips, geometry, counts, means, sds)
        values[site] = geometry.scaled ? value : value / 2
        compute_var && (variances[site] = variance_by_richness[length(tips)])
    end
    return PSVResult(psv = values, richness = richness, variance = variances, scale_vcv = geometry.scaled)
end

function psr(community::AbstractMatrix{<:Real}, tree::CompactTree; kwargs...)
    result = psv(community, tree; kwargs...)
    return (psr = result.psv .* result.richness, richness = result.richness,
        variance = result.variance .* result.richness .^ 2)
end

# Prepare weighted rows once for PSE/PSD, retaining CSC storage and native
# abundance types until normalization. Presence-only metrics need no weights.
function _species_abundance_rows(community, tree, tip_labels)
    _validate_community_values(community)
    size(community, 2) == tree.ntips ||
        throw(ArgumentError("community matrix must have one column per tree tip"))
    perm = _tools_tip_permutation(tree, tip_labels)
    n, T = size(community, 1), eltype(community)
    tips, raw = [Int32[] for _ in 1:n], [T[] for _ in 1:n]
    _community_each_positive(community) do row, col, value
        push!(tips[row], Int32(perm[col]))
        push!(raw[row], value)
    end
    weights = [Float64[] for _ in 1:n]
    precise = Union{Nothing, Vector{BigFloat}}[nothing for _ in 1:n]
    uniform = falses(n)
    for row in 1:n
        isempty(raw[row]) && continue
        scale = maximum(raw[row])
        uniform[row] = all(==(scale), raw[row])
        weights[row] = [_community_distance_ratio(x, scale) for x in raw[row]]
        weights[row] ./= sum(weights[row])
        if any(weight -> weight < floatmin(Float64), weights[row])
            # Zero/subnormal weights lose relative precision in the pair
            # product, even when the final scaled result is representable.
            # Long branches or many tiny records can also recover a value.
            values = BigFloat.(raw[row])
            values ./= maximum(values)
            precise[row] = values ./ sum(values)
        end
    end
    return (; tips, weights, precise, uniform)
end

# With normalized weights w, PSE = 2*sum(wi*wj*Yij)/(1-1/richness).
# Raw geometry stores patristic rather than half-patristic contrasts.
# Nonnegative cross-subtree contributions replace trace minus quadratic form.
function _pse_moments!(tree, tips, weights, geometry, mass, means)
    fill!(mass, 0)
    fill!(means, 0)
    for (tip, weight) in zip(tips, weights)
        mass[tree.tip_ids[tip]] = weight
    end
    pairs = zero(eltype(mass))
    for node in tree.postorder
        a, m = mass[node], zero(eltype(means))
        for child in tree.children[node]
            b = mass[child]
            iszero(b) && continue
            mc = geometry.offsets[child] + geometry.factors[child] * means[child]
            if iszero(a)
                a, m = b, mc
            else
                cross = geometry.scaled ? m + (1 - m) * mc : m + mc
                # Multiply the contrast before the two small weights.
                pairs += (cross * a) * b
                total = a + b
                m = (a / total) * m + (b / total) * mc
                a = total
            end
        end
        mass[node], means[node] = a, m
    end
    return (geometry.scaled ? 2pairs : pairs) / (1 - 1 / length(tips))
end

function _pse_values(rows, tree, geometry, mass, means; variability = nothing)
    result = fill(NaN, length(rows.tips))
    precise_mass, precise_means = BigFloat[], BigFloat[]
    for (site, tips) in enumerate(rows.tips)
        length(tips) > 1 || continue
        if variability !== nothing && rows.uniform[site]
            result[site] = variability[site]
        elseif rows.precise[site] !== nothing
            if isempty(precise_mass)
                precise_mass, precise_means = zeros(BigFloat, tree.nnodes), zeros(BigFloat, tree.nnodes)
            end
            result[site] = Float64(_pse_moments!(tree, tips, rows.precise[site],
                geometry, precise_mass, precise_means))
        else
            result[site] = _pse_moments!(tree, tips, rows.weights[site], geometry, mass, means)
        end
    end
    return result
end

"""
    pse(community, tree; scale_vcv=true, tip_labels=nothing)

Phylogenetic species evenness using normalized abundances and Brownian
correlation (default), or raw covariance. Empty/singleton sites give NaN.
Finite abundances are scaled before summation. Branch lengths must be finite,
nonnegative and have a finite sum; scaled calls need positive tip heights.

Uses sparse records and weighted subtree pair contrasts, without VCV or site
submatrices. Sparse calls reuse an observed connecting tree and can extract
smaller site trees; contracted affine paths preserve the shared root stem.
With equal positive abundances PSE equals PSV. Inputs are unchanged.
"""
function pse(community::AbstractMatrix{<:Real}, tree::CompactTree;
    scale_vcv::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    _validate_community_branch_lengths(tree)
    rows = _species_abundance_rows(community, tree, tip_labels)
    geometry = _psv_geometry(tree, scale_vcv)
    plan = _species_plan(rows.tips, tree, geometry)
    values = _species_pse_values(rows, plan)
    return (pse = values, richness = length.(rows.tips))
end

# Nearest OTHER species under 1-Cij. In correlation mode an edge lifts
# deficits by offset + factor*deficit in either direction. In raw mode all
# descendant candidates share covariance equal to the parent's height, while
# outside candidates retain their covariance on the downward pass.
function _psc_value!(workspace, tree, tips, geometry)
    d1, d2 = workspace.first_distance, workspace.second_distance
    t1, t2 = workspace.first_tip, workspace.second_tip
    fill!(d1, Inf)
    fill!(d2, Inf)
    fill!(t1, 0)
    fill!(t2, 0)
    for tip in tips
        node = tree.tip_ids[tip]
        d1[node] = geometry.scaled ? 0.0 : 1.0 - geometry.heights[node]
        t1[node] = node
    end
    for node in tree.postorder
        parent = tree.parent_of_node[node]
        parent == 0 && continue
        first = geometry.scaled ? geometry.offsets[node] + geometry.factors[node] * d1[node] :
            1.0 - geometry.heights[parent]
        second = geometry.scaled ? geometry.offsets[node] + geometry.factors[node] * d2[node] :
            1.0 - geometry.heights[parent]
        _comdistnt_insert!(workspace, parent, first, t1[node])
        _comdistnt_insert!(workspace, parent, second, t2[node])
    end
    for node in tree.preorder
        parent = tree.parent_of_node[node]
        parent == 0 && continue
        first = geometry.scaled ? geometry.down_offsets[node] + geometry.factors[node] * d1[parent] : d1[parent]
        second = geometry.scaled ? geometry.down_offsets[node] + geometry.factors[node] * d2[parent] : d2[parent]
        _comdistnt_insert!(workspace, node, first, t1[parent])
        _comdistnt_insert!(workspace, node, second, t2[parent])
    end
    n, value = 0, 0.0
    for tip in tips
        node = tree.tip_ids[tip]
        nearest = t1[node] == node ? d2[node] : d1[node]
        n, value = _psv_merge_mean(n, value, 1, nearest)
    end
    return value
end

function _psc_values(rows, tree, geometry)
    result = fill(NaN, length(rows))
    workspace = _ComdistntWorkspace(tree.nnodes)
    for (site, tips) in enumerate(rows)
        length(tips) > 1 || continue
        result[site] = _psc_value!(workspace, tree, tips, geometry)
    end
    return result
end

"""
    psc(community, tree; scale_vcv=true, tip_labels=nothing)

Phylogenetic species clustering: mean over present tips of
`1 - maximum(Cij for j != i)`. C is Brownian correlation by default, or raw
covariance when `scale_vcv=false`; raw PSC may be negative and depends on
branch-length units. Positive values denote presence. Empty/singleton sites
give NaN, and zero-length branches can give zero contrast to another species.

Two passes on full or observed connecting trees retain the nearest two
distinct identities without VCV or site submatrices. Original heights and
stable affine paths preserve raw/scaled definitions under contraction.
Branch lengths must be finite/nonnegative with finite sum; scaled calls
need positive tip heights. Inputs are not mutated.
"""
function psc(community::AbstractMatrix{<:Real}, tree::CompactTree;
    scale_vcv::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    _validate_community_branch_lengths(tree)
    rows = _community_tip_rows(community, tree, tip_labels)
    geometry = _psv_geometry(tree, scale_vcv; need_heights = !scale_vcv)
    plan = _species_plan(rows, tree, geometry)
    return (psc = _species_psc_values(plan), richness = length.(rows))
end

"""
    psd(community, tree; compute_var=true, scale_vcv=true, tip_labels=nothing)

Combined PSV, variance, PSC, PSR, PSR variance, PSE and richness, in input
site order. Shares label/record preparation, tree geometry, observed
connecting trees and workspace; derives PSR from PSV. Variance uses the
complete original species pool, including unobserved tips. No full VCV or
site submatrix is built, including with `compute_var=true`.

Keeps the single-metric definitions, validation and NaN rules. Inputs are
unchanged; each output field owns its array. Equal-abundance PSE reuses PSV.
"""
function psd(community::AbstractMatrix{<:Real}, tree;
    compute_var::Bool = true,
    scale_vcv::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    _validate_community_branch_lengths(tree)
    rows = _species_abundance_rows(community, tree, tip_labels)
    geometry = _psv_geometry(tree, scale_vcv; need_heights = !scale_vcv)
    plan = _species_plan(rows.tips, tree, geometry)
    if plan.index === nothing && plan.tree isa CompactTree
        # The original sequential route reuses its full-pool variance scratch
        # for site moments; allocating a second full workspace would regress
        # dense inputs without reducing their tree traversal.
        counts, means, scratch = zeros(Int, tree.nnodes), zeros(tree.nnodes), zeros(tree.nnodes)
        v = _psv_values(rows.tips, tree, geometry, counts, means, scratch; compute_var)
        evenness = _pse_values(rows, tree, geometry, scratch, means; variability = v.psv)
        clustering = _psc_values(rows.tips, tree, geometry)
        return (psv = v.psv, var_psv = v.variance, psc = clustering,
            psr = v.psv .* v.richness, var_psr = v.variance .* v.richness .^ 2,
            pse = evenness, richness = v.richness)
    end
    variance_by_richness = compute_var ? _psv_variances(tree, geometry,
        zeros(Int, tree.nnodes), zeros(tree.nnodes), zeros(tree.nnodes)) : Float64[]
    return _species_psd_values(rows, plan, variance_by_richness)
end
