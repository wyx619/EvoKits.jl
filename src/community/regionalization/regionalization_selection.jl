"""
    LinkageSelectionResult

`methods` and Pearson `correlations` between original and cophenetic distances,
in requested order. `best_method` is the first finite maximum, or `nothing` if
all correlations are undefined. `component` records the source beta component.
"""
struct LinkageSelectionResult
    methods::Vector{Symbol}
    correlations::Vector{Float64}
    best_method::Union{Nothing, Symbol}
    component::Symbol
end

"""
    OptimalPhyloregionResult

`df` contains `k`, explained fraction `ev`, within sum `totwss`, between sum
`totbss`, and total sum `tss`. `optimal_k` and selected `inc_threshold` /
`ev_threshold` are `nothing` when no recommendation is available; `reason` is
`:selected`, `:zero_dispersion`, or `:thresholds_not_met`. Method, component and
the internal normalization `distance_scale` are recorded. Sums use original
squared-distance units and may overflow to Inf or underflow to zero; `ev` and
threshold decisions are computed on scaled distances to avoid that overflow.
"""
struct OptimalPhyloregionResult
    df::DataFrame
    optimal_k::Union{Nothing, Int}
    inc_threshold::Union{Nothing, Float64}
    ev_threshold::Union{Nothing, Float64}
    reason::Symbol
    distance_scale::Float64
    method::Symbol
    component::Symbol
end

@inline function _region_compensated_add(value, correction, increment)
    adjusted = increment - correction
    updated = value + adjusted
    return updated, (updated - value) - adjusted
end

# Each original unordered pair is visited once, at its least common ancestor.
# Linked lists represent active leaf sets in O(n) storage, including combs.
# No cophenetic matrix or copies of descendant vectors are constructed.
function _region_tree_statistics(packed, merge, heights, n)
    scale = isempty(packed) ? 0.0 : maximum(packed)
    height_scale = isempty(heights) ? 0.0 : maximum(heights)
    heads, tails = collect(1:(2n - 1)), collect(1:(2n - 1))
    following = zeros(Int, n)
    sizes, squares = ones(Int, 2n - 1), zeros(2n - 1)
    within = zeros(n)
    count = 0
    mean_x, mean_y, variance_x, variance_y, covariance = 0.0, 0.0, 0.0, 0.0, 0.0
    total_within, correction = 0.0, 0.0
    for step in 1:(n - 1)
        a, b = merge[step, 1], merge[step, 2]
        left, right = a < 0 ? -a : n + a, b < 0 ? -b : n + b
        y = height_scale == 0 ? 0.0 : heights[step] / height_scale
        cross, cross_correction = 0.0, 0.0
        i = heads[left]
        while i != 0
            j = heads[right]
            while j != 0
                x = scale == 0 ? 0.0 : packed[_region_pair_index(i, j)] / scale
                cross, cross_correction = _region_compensated_add(cross, cross_correction, x^2)
                count += 1
                dx, dy = x - mean_x, y - mean_y
                mean_x += dx / count
                mean_y += dy / count
                variance_x += dx * (x - mean_x)
                variance_y += dy * (y - mean_y)
                covariance += dx * (y - mean_y)
                j = following[j]
            end
            i = following[i]
        end
        node = n + step
        sizes[node] = sizes[left] + sizes[right]
        squares[node] = squares[left] + squares[right] + cross
        change = squares[node] / sizes[node] - squares[left] / sizes[left] - squares[right] / sizes[right]
        total_within, correction = _region_compensated_add(total_within, correction, change)
        within[n - step] = max(0.0, total_within)
        heads[node], tails[node] = heads[left], tails[right]
        following[tails[left]] = heads[right]
    end
    total = squares[end] / n
    within[1] = total # exact k=1 boundary, without accumulation roundoff.
    correlation = variance_x > 0 && variance_y > 0 ?
        clamp(covariance / sqrt(variance_x * variance_y), -1.0, 1.0) : NaN
    return correlation, within, total, scale
end

function _region_method_list(methods)
    methods isa Symbol && throw(ArgumentError("methods must be a collection of method symbols"))
    values = Symbol[]
    for method in methods
        method isa Symbol || throw(ArgumentError("methods must contain symbols"))
        push!(values, _region_method(method))
    end
    isempty(values) && throw(ArgumentError("methods must not be empty"))
    allunique(values) || throw(ArgumentError("methods must be unique after alias normalization"))
    return values
end

function _select_linkage(packed, n, component, methods, context = nothing)
    values = _region_method_list(methods)
    correlations = Float64[]
    best, best_value = nothing, -Inf
    for method in values
        merge, heights = context !== nothing && method === context.method ?
            (context.merge, context.heights) : _region_hierarchy!(copy(packed), n, method)
        correlation, _, _, _ = _region_tree_statistics(packed, merge, heights, n)
        push!(correlations, correlation)
        if isfinite(correlation) && correlation > best_value
            best, best_value = method, correlation
        end
    end
    return LinkageSelectionResult(values, correlations, best, component)
end

"""
    select_linkage(distances; methods=(:average, :single, :complete,
        :ward_D, :ward_D2, :mcquitty, :median, :centroid))
    select_linkage(beta; component=:total, methods=...)
    select_linkage(context::PhyloregionContext; methods=...)

Compare original distances with each hierarchy's cophenetic distances using
Pearson correlation. Larger values indicate better distance preservation, not
ecological correctness or an optimal k. Tied maxima choose the first requested
method. Constant distances/heights or fewer than two pairs give NaN; if every
score is undefined, `best_method` is `nothing`. Context input reuses its packed
distances and the existing tree for its method. Methods run serially; no full
cophenetic matrix is allocated. Ward's variance interpretation requires
Euclidean distances; centroid/median use R's recurrences on the supplied values.
"""
function select_linkage(distances::AbstractMatrix{<:Real}; methods = _REGION_METHODS)
    packed = _region_distance_buffer(distances)
    return _select_linkage(packed, size(distances, 1), :distance, methods)
end

function select_linkage(beta::Union{BetaDissResult, PhylobetaResult};
    component::Symbol = :total, methods = _REGION_METHODS,
)
    distances = _region_beta_distances(beta, component)
    packed = _region_distance_buffer(distances)
    return _select_linkage(packed, size(distances, 1), component, methods)
end

function select_linkage(context::PhyloregionContext; methods = _REGION_METHODS)
    return _select_linkage(context.distances, length(context.sample_labels),
        context.component, methods, context)
end

function _region_thresholds(values, name; rev = false)
    thresholds = Float64[]
    for value in values
        value isa Real && !(value isa Bool) && isfinite(value) && 0 <= value <= 1 ||
            throw(ArgumentError("$name must contain finite numbers between 0 and 1"))
        push!(thresholds, Float64(value))
    end
    isempty(thresholds) && throw(ArgumentError("$name must not be empty"))
    return sort!(unique(thresholds); rev)
end

function _region_original_ss(value, scale)
    value == 0 && return 0.0
    root = sqrt(abs(value)) * scale
    return copysign(root * root, value)
end

function _region_optimal(context, k, inc_thresholds, ev_thresholds, precision)
    n = length(context.sample_labels)
    upper = k === nothing ? min(n, 20) : _region_validate_k(k, n)
    precision isa Integer && !(precision isa Bool) && 0 <= precision <= 15 ||
        throw(ArgumentError("precision must be an integer between 0 and 15"))
    increments = _region_thresholds(inc_thresholds, "inc_thresholds")
    fractions = _region_thresholds(ev_thresholds, "ev_thresholds"; rev = true)
    _, within, total, scale = _region_tree_statistics(context.distances, context.merge, context.heights, n)
    ev = total > 0 ? [1 - within[i] / total for i in 1:upper] : fill(NaN, upper)
    table = DataFrame(k = collect(1:upper), ev = ev,
        totwss = [_region_original_ss(within[i], scale) for i in 1:upper],
        totbss = [_region_original_ss(total - within[i], scale) for i in 1:upper],
        tss = fill(_region_original_ss(total, scale), upper))
    selected, chosen_inc, chosen_ev = nothing, nothing, nothing
    reason = total == 0 ? :zero_dispersion : :thresholds_not_met
    tolerance = 10.0^(-precision)
    # Preserve the saved phyloregion elbow.batch threshold order and tolerance.
    for ev_threshold in fractions, inc_threshold in increments
        first_ev = findfirst(value -> isfinite(value) && value >= ev_threshold - tolerance, ev)
        first_inc = findfirst(i -> isfinite(ev[i + 1] - ev[i]) &&
            ev[i + 1] - ev[i] <= inc_threshold + tolerance, 1:(upper - 1))
        if first_ev !== nothing && first_inc !== nothing
            selected = max(first_ev, first_inc)
            chosen_inc, chosen_ev = inc_threshold, ev_threshold
            reason = :selected
            break
        end
    end
    return OptimalPhyloregionResult(table, selected, chosen_inc, chosen_ev, reason,
        scale, context.method, context.component)
end

"""
    optimal_phyloregion(context; k=nothing,
        inc_thresholds=(0.01, 0.05, 0.1),
        ev_thresholds=(0.95, 0.9, 0.8, 0.75, 0.67, 0.5, 0.33), precision=3)
    optimal_phyloregion(distances; method=:average, sample_labels=nothing, ...)
    optimal_phyloregion(beta; component=:total, method=:average, ...)

Evaluate cuts 1:k (default min(n,20)) on one hierarchy. For a cluster C, WSS is
sum(d[i,j]^2 for unordered pairs in C)/length(C); TSS uses all sites, and
`ev = (TSS - total_WSS)/TSS`. For non-Euclidean dissimilarities these are
distance-based scores: ev need not be monotone or lie in [0,1].

Use saved phyloregion threshold semantics: try EV thresholds in descending
order, then increment thresholds ascending. For each combination, find the
first k meeting the EV threshold and the first k whose next increment meets
the increment threshold, with tolerance 10^(-precision); recommend their max.
The last evaluated k has no next increment. This is a threshold heuristic,
not a curvature optimizer or a guarantee of a biologically correct partition.
Zero dispersion or unmet thresholds return `optimal_k=nothing` plus the full
table and reason. Squared-distance sums may overflow/underflow in original
units, but ev/decisions use normalized distances. All pairs are visited once
using O(n) auxiliary storage; a context avoids rebuilding the hierarchy.
"""
function optimal_phyloregion(context::PhyloregionContext; k = nothing,
    inc_thresholds = (0.01, 0.05, 0.1),
    ev_thresholds = (0.95, 0.9, 0.8, 0.75, 0.67, 0.5, 0.33), precision = 3,
)
    return _region_optimal(context, k, inc_thresholds, ev_thresholds, precision)
end

function optimal_phyloregion(distances::AbstractMatrix{<:Real};
    method::Symbol = :average, sample_labels = nothing, kwargs...,
)
    return optimal_phyloregion(prepare_phyloregion(distances; method, sample_labels); kwargs...)
end

function optimal_phyloregion(beta::Union{BetaDissResult, PhylobetaResult};
    method::Symbol = :average, component::Symbol = :total, kwargs...,
)
    return optimal_phyloregion(prepare_phyloregion(beta; method, component); kwargs...)
end
