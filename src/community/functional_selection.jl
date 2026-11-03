"""Candidate functional-group fits and an optional reusable elbow-selected grouping."""
struct FunctionalGroupSelectionResult
    k::Vector{Int}
    objectives::Vector{Float64}
    iterations::Vector{Int}
    converged::BitVector
    optimal_k::Union{Nothing, Int}
    grouping::Union{Nothing, FunctionalGroupingResult}
    reason::Symbol
    low::Float64
    max_fraction::Float64
end

# The lbow heuristic from phyloregion, applied to one comparable objective curve.
# Scaling prevents overflow in the normalization; degenerate curves are explicit.
function _functional_elbow(objectives; low=0.08, max_fraction=0.9)
    largest = maximum(objectives)
    largest > 0 || return nothing, :zero_objective
    weights = objectives ./ largest
    weights ./= sum(weights)
    threshold = Float64(low)
    while threshold >= maximum(weights)
        threshold /= 2
    end
    small = findall(x -> threshold / 8 < x < threshold, weights)
    index = if isempty(small)
        length(weights)
    elseif length(small) == 1
        only(small)
    else
        drops = [abs(weights[small[i + 1]] - weights[small[i]]) / weights[small[i]]
            for i in 1:(length(small) - 1)]
        # R uses the number of all above-threshold entries as its offset.
        argmax(drops) + count(>=(threshold), weights)
    end
    if sum(@view weights[1:index]) > max_fraction
        index = findfirst(>(max_fraction), cumsum(weights)) - 1
    end
    index >= 1 || return nothing, :thresholds_not_met
    return index, :ok
end

function _optimal_functional_groups(table, selected; k=2:10, low=0.08, max_fraction=0.9,
    species_col=:species, trait_cols=nothing, categorical_traits=(), missing_traits=:error,
    numeric_scale=:sd, trait_weights=nothing, gamma=:auto, nstart=10, maxiter=100,
    rng=Random.default_rng(),
)
    (k isa AbstractVector || k isa AbstractRange || k isa Tuple) ||
        throw(ArgumentError("k must be a collection of candidate group counts"))
    all(x -> x isa Integer && !(x isa Bool) && 1 <= x <= typemax(Int), k) ||
        throw(ArgumentError("candidate k values must be positive Int-sized integers"))
    counts = sort!(Int.(collect(k)))
    length(counts) >= 2 && allunique(counts) ||
        throw(ArgumentError("at least two distinct candidate k values are required"))
    for (value, name) in ((nstart, "nstart"), (maxiter, "maxiter"))
        value isa Integer && !(value isa Bool) && 1 <= value <= typemax(Int) ||
            throw(ArgumentError("$name must be a positive Int-sized integer"))
    end
    for (value, name) in ((low, "low"), (max_fraction, "max_fraction"))
        value isa Real && !(value isa Bool) && isfinite(value) && 0 < value < 1 &&
            0 < Float64(value) < 1 || throw(ArgumentError("$name must be between zero and one"))
    end
    work = _functional_workspace(table, selected; species_col, trait_cols, categorical_traits,
        missing_traits, numeric_scale, trait_weights, gamma)
    last(counts) <= length(work.species) ||
        throw(ArgumentError("candidate k exceeds the complete-trait species count"))
    fits = [_functional_fit_prepared(work; k=count, nstart, maxiter, rng) for count in counts]
    objectives = [fit.objective for fit in fits]
    converged = BitVector(fit.converged for fit in fits)
    index, reason = if !all(converged)
        nothing, :nonconverged_candidates
    elseif any(objectives[i + 1] > objectives[i] + 64eps(Float64) * max(objectives[i], objectives[i + 1])
        for i in 1:(length(objectives) - 1))
        nothing, :nonmonotone_objective
    else
        _functional_elbow(objectives; low, max_fraction)
    end
    return FunctionalGroupSelectionResult(counts, objectives, [fit.iterations for fit in fits],
        converged, index === nothing ? nothing : counts[index],
        index === nothing ? nothing : fits[index], reason, Float64(low), Float64(max_fraction))
end

"""
    optimal_functional_groups(traits; k=2:10, low=0.08, max_fraction=0.9, kwargs...)
    optimal_functional_groups(data::CommunityData, traits; kwargs...)

Fit candidate k-prototypes group counts in ascending order, preparing labels,
scales, weights and gamma once. Other keywords follow `functional_groups`.
The supplied RNG is consumed sequentially across candidates and starts.
Returns all objectives and convergence diagnostics, `optimal_k`, `grouping`,
and `reason`. Pass a selected grouping to `collapse_traits` or `functional_beta`.

Selection follows phyloregion's lbow heuristic, not a statistical optimum.
Zero objectives, nonconverged fits, increasing objective curves, or a failed
threshold return `nothing` for the recommendation and grouping. Candidate
counts must be distinct, include at least two values and not exceed the number
of complete-trait species. Candidate groupings are temporarily retained until
selection; the result retains only the selected one.
"""
optimal_functional_groups(traits::AbstractDataFrame; kwargs...) =
    _optimal_functional_groups(traits, nothing; kwargs...)

function optimal_functional_groups(data::CommunityData, traits::AbstractDataFrame; kwargs...)
    _validate_community_data(data)
    return _optimal_functional_groups(traits, data.species_labels; kwargs...)
end
