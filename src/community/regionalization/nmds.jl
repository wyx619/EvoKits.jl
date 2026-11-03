"""
    NMDSResult

Global nonmetric ordination returned by [`nmds`](@ref). Rows of `points` follow
`sample_labels`; scores are centered, rotated to principal axes and have unit
RMS distance to the centroid. `stress` is Kruskal stress-1, on the 0:1 scale.

`dissimilarities`, `distances` and `disparities` contain original input,
Euclidean score distances and monotone fitted distances, respectively, packed
in upper-triangle column order: (1,2), (1,3), (2,3), (1,4), ... . The latter
two use the score scale. `stress_history` includes the initial stress and each
accepted iteration of the winning start. `iterations`, `converged`, `reason`
refer to that start; `start_stress`, `start_iterations`, `start_converged` and
`start_reason` describe all starts, including the winner `best_start`.
Convergence reports a numerical stopping rule, not a proven global optimum.
Arrays are independently owned. `component` records the input beta component,
or :distance; `index_family` is nothing for plain or regional distances.
"""
struct NMDSResult
    sample_labels::Vector{String}
    points::Matrix{Float64}
    stress::Float64
    dissimilarities::Vector{Float64}
    distances::Vector{Float64}
    disparities::Vector{Float64}
    iterations::Int
    converged::Bool
    reason::Symbol
    best_start::Int
    start_stress::Vector{Float64}
    start_iterations::Vector{Int}
    start_converged::Vector{Bool}
    start_reason::Vector{Symbol}
    stress_history::Vector{Float64}
    dims::Int
    ties::Symbol
    component::Symbol
    index_family::Union{Nothing, Symbol}
end

struct _NMDSWorkspace
    order::Vector{Int}
    groups::Vector{UnitRange{Int}}
    means::Vector{Float64}
    weights::Vector{Int}
    ends::Vector{Int}
    distances::Vector{Float64}
    fitted::Vector{Float64}
end

function _NMDSWorkspace(dissimilarities)
    order = sortperm(dissimilarities)
    groups = UnitRange{Int}[]
    first = 1
    for last in eachindex(order)
        if last == length(order) || dissimilarities[order[last]] != dissimilarities[order[last + 1]]
            push!(groups, first:last)
            first = last + 1
        end
    end
    m = length(order)
    return _NMDSWorkspace(order, groups, zeros(m), zeros(Int, m), zeros(Int, m), zeros(m), zeros(m))
end

# Weighted pool-adjacent-violators regression. Weak ties sort current score
# distances within each exact input tie; strong ties start as indivisible
# weighted blocks. All scratch arrays are reused between iterations/starts.
function _nmds_disparities!(workspace, ties)
    (; order, groups, means, weights, ends, distances, fitted) = workspace
    if ties === :weak
        for group in groups
            length(group) > 1 && sort!(view(order, group); alg = QuickSort, by = i -> (distances[i], i))
        end
    end
    blocks = 0
    for group in (ties === :weak ? (i:i for i in eachindex(order)) : groups)
        weight = length(group)
        value = 0.0
        for i in group
            value += distances[order[i]] / weight
        end
        blocks += 1
        means[blocks], weights[blocks], ends[blocks] = value, weight, last(group)
        while blocks > 1 && means[blocks - 1] > means[blocks]
            a, b = weights[blocks - 1], weights[blocks]
            means[blocks - 1] += (means[blocks] - means[blocks - 1]) * (b / (a + b))
            weights[blocks - 1] = a + b
            ends[blocks - 1] = ends[blocks]
            blocks -= 1
        end
    end
    first = 1
    for block in 1:blocks
        for i in first:ends[block]
            fitted[order[i]] = means[block]
        end
        first = ends[block] + 1
    end
    return fitted
end

function _nmds_distances!(dest, points)
    n, dims = size(points)
    pair = 0
    for j in 2:n, i in 1:(j - 1)
        value = 0.0
        for axis in 1:dims
            value = hypot(value, points[i, axis] - points[j, axis])
        end
        pair += 1
        dest[pair] = value
    end
    return dest
end

function _nmds_stress!(workspace, points, ties)
    distances = _nmds_distances!(workspace.distances, points)
    fitted = _nmds_disparities!(workspace, ties)
    numerator, denominator = 0.0, 0.0
    for i in eachindex(distances)
        numerator += abs2(distances[i] - fitted[i])
        denominator += abs2(distances[i])
    end
    return denominator > 0 ? sqrt(numerator / denominator) : NaN
end

# Rescale before centering/squaring so finite extreme initial units cannot
# overflow or underflow. A completely collapsed configuration is inadmissible.
function _nmds_normalize!(points)
    scale = maximum(abs, points)
    scale > 0 || return false
    points ./= scale
    n, dims = size(points)
    for axis in 1:dims
        center = sum(view(points, :, axis)) / n
        for i in 1:n
            points[i, axis] -= center
        end
    end
    rms = norm(points) / sqrt(n)
    rms > 0 || return false
    points ./= rms
    return true
end

function _nmds_initial_matrix(initial, n, dims)
    initial isa AbstractMatrix{<:Real} && size(initial) == (n, dims) ||
        throw(ArgumentError("init must be :pcoa, :random or an nsites-by-dims real matrix"))
    Base.require_one_based_indexing(initial)
    all(isfinite, initial) || throw(ArgumentError("initial coordinates must be finite"))
    points = Matrix{Float64}(initial)
    all(isfinite, points) && all(pair -> pair[1] == 0 || pair[2] != 0, zip(initial, points)) ||
        throw(ArgumentError("initial coordinates must be representable as finite Float64 without underflow to zero"))
    _nmds_normalize!(points) || throw(ArgumentError("initial coordinates must not all coincide"))
    return points
end

function _nmds_pcoa(dissimilarities, n, dims)
    scale = maximum(dissimilarities)
    gram = zeros(n, n)
    for j in 2:n, i in 1:(j - 1)
        value = abs2(dissimilarities[_region_pair_index(i, j)] / scale)
        gram[i, j] = gram[j, i] = value
    end
    row_means = vec(sum(gram; dims = 2)) ./ n
    grand_mean = sum(row_means) / n
    for j in 1:n, i in 1:n
        gram[i, j] = -0.5 * (gram[i, j] - row_means[i] - row_means[j] + grand_mean)
    end
    decomposition = eigen!(Symmetric(gram))
    points = zeros(n, dims)
    for axis in 1:dims
        index = n - axis + 1
        value = decomposition.values[index]
        if value > 0
            for i in 1:n
                points[i, axis] = decomposition.vectors[i, index] * sqrt(value)
            end
        end
    end
    _nmds_normalize!(points) || throw(ArgumentError("PCoA initialization collapsed; supply init=:random or explicit coordinates"))
    return points
end

function _nmds_random!(points, rng)
    rand!(rng, points)
    points .-= 0.5
    _nmds_normalize!(points) || throw(ArgumentError("random initialization collapsed; supply a nonconstant RNG or explicit coordinates"))
    return points
end

# Complete equal-weight Guttman transform followed by unit-RMS normalization.
# Fixed score norm makes the stress-1 denominator constant; normalization of
# B(X)*X maximizes its linear majorization term on that sphere. Compute
# differences/distance before multiplying by the disparity to avoid overflow
# for nearly coincident points. A coincident pair contributes zero.
function _nmds_step!(dest, points, workspace)
    fill!(dest, 0.0)
    n, dims = size(points)
    pair = 0
    for j in 2:n, i in 1:(j - 1)
        pair += 1
        distance = workspace.distances[pair]
        distance > 0 || continue
        fitted = workspace.fitted[pair]
        for axis in 1:dims
            value = ((points[i, axis] - points[j, axis]) / distance) * fitted
            dest[i, axis] += value
            dest[j, axis] -= value
        end
    end
    return _nmds_normalize!(dest)
end

function _nmds_run!(points, candidate, workspace, ties, maxiter, tol, stress_tol)
    stress = _nmds_stress!(workspace, points, ties)
    history = [stress]
    stress <= stress_tol && return (; points, stress, iterations = 0, converged = true, reason = :stress_tol, history)
    for iteration in 1:maxiter
        _nmds_step!(candidate, points, workspace) || return (
            ; points, stress, iterations = iteration - 1, converged = false, reason = :collapsed, history)
        proposed = _nmds_stress!(workspace, candidate, ties)
        if !isfinite(proposed) || proposed > stress
            # Keep the last accepted scores and restore its pair work arrays.
            _nmds_stress!(workspace, points, ties)
            small = isfinite(proposed) && proposed - stress <= 64eps(Float64)
            return (; points, stress, iterations = iteration - 1,
                converged = small, reason = small ? :numerical_stall : :stress_increase, history)
        end
        previous = stress
        stress = proposed
        points, candidate = candidate, points
        push!(history, stress)
        reason = stress <= stress_tol ? :stress_tol :
            previous - stress <= tol * max(previous, eps(Float64)) ? :relative_tol : :continue
        reason === :continue || return (; points, stress, iterations = iteration, converged = true, reason, history)
    end
    return (; points, stress, iterations = maxiter, converged = false, reason = :maxiter, history)
end

function _nmds_integer(value, name, minimum, maximum = typemax(Int))
    value isa Integer && !(value isa Bool) && minimum <= value <= maximum ||
        throw(ArgumentError("$name must be an integer between $minimum and $maximum"))
    return Int(value)
end

function _nmds_tolerance(value, name)
    value isa Real && !(value isa Bool) && isfinite(value) && 0 <= value < 1 ||
        throw(ArgumentError("$name must be finite and between zero (inclusive) and one (exclusive)"))
    return Float64(value)
end

function _nmds(packed, labels, component, family;
    dims = 2, nstart = 20, maxiter = 1000, init = :pcoa, ties::Symbol = :weak,
    tol = 1e-7, stress_tol = 1e-4, rng::AbstractRNG = Random.default_rng(),
)
    n = length(labels)
    n > 1 || throw(ArgumentError("NMDS requires at least two sites"))
    length(packed) == n * (n - 1) ÷ 2 && all(x -> isfinite(x) && x >= 0, packed) ||
        throw(ArgumentError("NMDS requires complete finite nonnegative distances matching the labels"))
    k = _nmds_integer(dims, "dims", 1, n - 1)
    starts = _nmds_integer(nstart, "nstart", 1)
    limit = _nmds_integer(maxiter, "maxiter", 0)
    tolerance, minimum_stress = _nmds_tolerance(tol, "tol"), _nmds_tolerance(stress_tol, "stress_tol")
    ties in (:weak, :strong) || throw(ArgumentError("ties must be :weak or :strong"))
    maximum(packed) > 0 || throw(ArgumentError("NMDS requires at least one positive dissimilarity"))
    first_points = init === :pcoa ? _nmds_pcoa(packed, n, k) :
        init === :random ? _nmds_random!(zeros(n, k), rng) : _nmds_initial_matrix(init, n, k)
    workspace = _NMDSWorkspace(packed)
    points, candidate = first_points, similar(first_points)
    scores, iterations = zeros(starts), zeros(Int, starts)
    converged, reasons = fill(false, starts), fill(:maxiter, starts)
    best_points, best_history = similar(first_points), Float64[]
    best_start, best_stress = 0, Inf
    for start in 1:starts
        start > 1 && _nmds_random!(points, rng)
        run = _nmds_run!(points, candidate, workspace, ties, limit, tolerance, minimum_stress)
        scores[start], iterations[start] = run.stress, run.iterations
        converged[start], reasons[start] = run.converged, run.reason
        if run.stress < best_stress
            best_start, best_stress = start, run.stress
            copyto!(best_points, run.points)
            best_history = run.history
        end
        # The solver swaps its two score buffers; retain both without aliasing.
        points, candidate = run.points === points ? (points, candidate) : (candidate, points)
    end
    # Rotate only the winner. Principal-axis signs use the largest-magnitude
    # score as a reproducible positive anchor; tied eigenaxes remain arbitrary.
    rotation = svd(best_points; full = false).V
    best_points = best_points * rotation
    for axis in 1:k
        anchor = argmax(abs.(view(best_points, :, axis)))
        best_points[anchor, axis] < 0 && (best_points[:, axis] .*= -1)
    end
    _nmds_stress!(workspace, best_points, ties)
    return NMDSResult(copy(labels), best_points, best_stress, copy(packed),
        copy(workspace.distances), copy(workspace.fitted), iterations[best_start],
        converged[best_start], reasons[best_start], best_start, scores, iterations,
        converged, reasons, best_history, k, ties, component, family)
end

"""
    nmds(distances; sample_labels=nothing, dims=2, nstart=20, maxiter=1000,
        init=:pcoa, ties=:weak, tol=1e-7, stress_tol=1e-4,
        rng=Random.default_rng())
    nmds(beta; component=:total, ...)
    nmds(average::MeanDistResult; ...)
    nmds(regions::PhyloregionResult; ...)
    nmds(context::PhyloregionContext; ...)

Global nonmetric multidimensional scaling with Kruskal stress-1:
`sqrt(sum((distance-disparity)^2)/sum(distance^2))`. Fits monotone disparities
with PAVA and uses normalized SMACOF majorization. Weak ties allow equal input
dissimilarities to have different fitted values; strong ties require identical
fitted values. Ties are exact, unlike vegan's absolute near-tie tolerance.

Input must be a complete, symmetric, finite, nonnegative distance matrix with
zero diagonal and at least one positive entry. NaN/Inf and selected-pair
results are rejected; no sites are silently removed or distances imputed.
Off-diagonal zeros participate in regression and do not force coincident
scores. Constant positive inputs under weak ties can have zero stress for any
configuration. Require 1 <= dims < nsites; NMDS generally needs many more sites
than dimensions for an informative fit.

`nstart` counts all starts. The first uses PCoA (`init=:pcoa`, default), random
scores (`:random`) or a supplied nsites-by-dims matrix. Other starts are random.
The same RNG controls all random starts; increasing nstart retains the earlier
random sequence. Each start runs independently, and the smallest final stress
wins (first on exact ties). PCoA uses positive eigenaxes, pads absent axes with
zeros, and requires O(nsites^3) work; :random avoids that decomposition.
Rank-deficient starts can remain in their initial subspace; random starts
explore additional directions.

Output scores are centered, unit RMS and rotated to principal axes. Positive
input units do not determine score units. Numerical convergence means stress
<= stress_tol, relative stress decrease <= tol, or a rounding-level stall.
`maxiter=0` evaluates starts without updating them. Inspect per-start diagnostics;
neither convergence nor multiple starts establishes a global optimum. Only the
best stress history is retained; all runs finish before selection. Iteration
storage is O(nsites^2 + nsites*dims), independent of nstart*maxiter; the history
adds O(maxiter) and diagnostics O(nstart).

Beta inputs select :total/:turnover/:nestedness and retain site labels.
`MeanDistResult` uses its selected component. `PhyloregionResult` ordains its
region_distances with labels "1", "2", ... in membership-ID order; a prepared
context ordains original site distances without rebuilding its clustering.
Functional beta inputs use their contained beta result. Inputs are unchanged.
This numeric core does not transform community data, fit species scores,
apply step-across/half-change scaling, or construct maps and colors.
"""
function nmds(distances::AbstractMatrix{<:Real}; sample_labels = nothing, kwargs...)
    labels = _validated_labels(sample_labels, size(distances, 1), "sample")
    return _nmds(_region_distance_buffer(distances), labels, :distance, nothing; kwargs...)
end

function nmds(beta::Union{BetaDissResult, PhylobetaResult}; component::Symbol = :total, kwargs...)
    _community_data_keywords(kwargs)
    _validate_community_beta_family(beta.index_family)
    packed = _region_distance_buffer(_region_beta_distances(beta, component))
    labels = _validated_labels(beta.sample_labels, length(beta.sample_labels), "sample")
    length(labels) * (length(labels) - 1) ÷ 2 == length(packed) ||
        throw(ArgumentError("beta matrix must match sample labels"))
    return _nmds(packed, labels, component, beta.index_family; kwargs...)
end

nmds(result::FunctionalBetaResult; kwargs...) = nmds(result.beta; kwargs...)

function nmds(result::MeanDistResult; kwargs...)
    _community_data_keywords(kwargs)
    result.pairs === nothing && result.distances isa AbstractMatrix ||
        throw(ArgumentError("NMDS requires a full mean distance matrix, not selected pairs"))
    n = length(result.sample_labels)
    size(result.distances) == (n, n) || throw(ArgumentError("mean distance matrix must match sample labels"))
    labels = _validated_labels(result.sample_labels, n, "sample")
    return _nmds(_region_distance_buffer(result.distances), labels, result.component, result.index_family; kwargs...)
end

function nmds(regions::PhyloregionResult; kwargs...)
    _community_data_keywords(kwargs)
    return _nmds(_region_distance_buffer(regions.region_distances), string.(1:regions.k),
        regions.component, nothing; kwargs...)
end

function nmds(context::PhyloregionContext; kwargs...)
    _community_data_keywords(kwargs)
    n = length(context.sample_labels)
    labels = _validated_labels(context.sample_labels, n, "sample")
    return _nmds(context.distances, labels, context.component, nothing; kwargs...)
end
