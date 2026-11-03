"""
    GoMResult

Fixed-k grade-of-membership fit returned by [`fitgom`](@ref). `omega` has sites
× k shape with row sums one; `theta` has species × k shape with column sums one.
Labels retain input order. Regions are ordered by decreasing total `omega`,
with stable ties. `membership` takes each site's first maximum; `cluster_sizes`
counts those hard assignments and can include zeroes. Mixture weights are not
classification confidence or posterior uncertainty.

`objective` is the unnormalized NEF penalized log likelihood, `loglikelihood`
omits the priors and multinomial constants, and `history` records monotonic
objectives for the chosen start, including initialization. `residual` is the
largest absolute change in either factor under one ordinary EM update.
`iterations`, `converged`, `reason` describe the chosen start; `start_*` fields
record all starts. Each iteration takes two EM steps and may accept a
safeguarded log-ratio extrapolation. All result arrays are independently owned.
"""
struct GoMResult
    sample_labels::Vector{String}
    species_labels::Vector{String}
    region_labels::Vector{String}
    omega::Matrix{Float64}
    theta::Matrix{Float64}
    k::Int
    membership::Vector{Int}
    cluster_sizes::Vector{Int}
    loglikelihood::Float64
    objective::Float64
    shape::Float64
    iterations::Int
    converged::Bool
    reason::Symbol
    residual::Float64
    history::Vector{Float64}
    nstart::Int
    maxiter::Int
    tol::Float64
    accelerate::Bool
    best_start::Int
    start_objectives::Vector{Float64}
    start_iterations::Vector{Int}
    start_converged::BitVector
    start_reasons::Vector{Symbol}
    start_residuals::Vector{Float64}
end

function _gom_integer(value, name)
    value isa Integer && !(value isa Bool) && 1 <= value <= typemax(Int) ||
        throw(ArgumentError("$name must be a positive Int-sized integer"))
    return Int(value)
end

"""
    fitgom(community; k, shape=nothing, nstart=5, maxiter=1000, tol=1e-8,
           accelerate=true, init_theta=nothing, init_omega=nothing,
           rng=Random.default_rng(), sample_labels=nothing, species_labels=nothing)
    fitgom(data::CommunityData; k, ...)

Fit a fixed-k grade-of-membership model, with `1 <= k <= min(sites, species)`.
For each site i and species j, q_ij = sum_r omega[i,r]*theta[j,r]. The scalar
`shape` defaults to 1/(k*species), matching maptpx's NEF objective:

    sum_ij x_ij*log(q_ij) + shape*sum(log, theta) + sum(log, omega)/k

`shape` is a strictly positive pseudocount, not shape minus one. Finite
nonnegative counts are used as supplied: binary data give presence-based fits;
integer abundances give multinomial fits. Fractional weights are allowed as a
weighted likelihood. Scaling counts changes their balance with the priors.
Empty sites are errors; all species columns, including all-zero columns, stay
in `theta` and in the default prior's species count. CSC input visits stored
entries only and skips explicit zeroes; no dense site × species or nnz × k
responsibility matrix is created. Work storage is O(nnz + (sites+species)*k).

Positive `init_theta` (species × k) and `init_omega` (sites × k) are normalized
and used for the first start. Unspecified factors and subsequent starts are
randomized using `rng`; the best objective wins, with the first start winning
exact ties. No input is mutated. Each iteration takes two EM steps and may
accept an improving extrapolation in log ratios; `accelerate=false` uses EM
alone. This optimizes the scalar-shape maptpx objective, using a different
optimizer from its conditional SQP/EM routine. Local optima can differ.

Convergence requires the maximum absolute one-step EM change to be at most
`tol`. Stopping reasons are `:tolerance`, `:maxiter`, or `:stalled` (no numerical
increase); the latter two need not be converged. No automatic k selection,
Bayes factors, geographic constraint, or uncertainty estimates are computed.
For lineage/functional models, pass an aggregated `.community`. Pass the
returned fit to [`indicators`](@ref) for KL indicator species.
"""
function fitgom(community::AbstractMatrix{<:Real}; k, shape=nothing, nstart=5,
    maxiter=1000, tol=1e-8, accelerate=true, init_theta=nothing, init_omega=nothing,
    rng::AbstractRNG=Random.default_rng(), sample_labels=nothing, species_labels=nothing)
    groups = _gom_integer(k, "k")
    starts = _gom_integer(nstart, "nstart")
    limit = _gom_integer(maxiter, "maxiter")
    groups <= min(size(community)...) ||
        throw(ArgumentError("k cannot exceed the number of sites or species"))
    tol isa Real && !(tol isa Bool) && isfinite(tol) && 0 < tol < 1 ||
        throw(ArgumentError("tol must be finite and between zero and one"))
    tolerance = Float64(tol)
    tolerance > 0 || throw(ArgumentError("tol underflows Float64"))
    accelerate isa Bool || throw(ArgumentError("accelerate must be a Bool"))
    samples = _validated_labels(sample_labels, size(community, 1), "sample")
    species = _validated_labels(species_labels, size(community, 2), "species")
    counts = _gom_counts(community)
    prior = shape === nothing ? 1 / (Float64(groups) * counts.p) : shape
    prior isa Real && !(prior isa Bool) && isfinite(prior) && prior > 0 ||
        throw(ArgumentError("shape must be a finite positive scalar"))
    alpha = Float64(prior)
    isfinite(alpha) && alpha > 0 && isfinite(counts.total + counts.p * alpha) &&
        alpha / (counts.total + counts.p * alpha) > 0 ||
        throw(ArgumentError("shape and counts cannot be combined in Float64 range; rescale counts or priors"))
    initial_theta = init_theta === nothing ? nothing :
        _gom_initial(init_theta, counts.p, groups, false, "init_theta")
    initial_omega = init_omega === nothing ? nothing :
        _gom_initial(init_omega, counts.n, groups, true, "init_omega")
    objectives, iterations, residuals = Float64[], Int[], Float64[]
    converged, reasons = BitVector(), Symbol[]
    best, best_start = nothing, 0
    for start in 1:starts
        theta, omega = if start == 1 && initial_theta !== nothing && initial_omega !== nothing
            copy(initial_theta), copy(initial_omega)
        else
            _gom_random_start(counts, groups, rng)
        end
        start == 1 && initial_theta !== nothing && (theta = copy(initial_theta))
        start == 1 && initial_omega !== nothing && (omega = copy(initial_omega))
        fit = _gom_fit(counts, theta, omega, alpha, limit, tolerance, accelerate)
        push!(objectives, fit.objective)
        push!(iterations, fit.iterations)
        push!(converged, fit.converged)
        push!(reasons, fit.reason)
        push!(residuals, fit.residual)
        if best === nothing || fit.objective > best.objective
            best, best_start = fit, start
        end
    end
    usage = vec(sum(best.omega; dims=1))
    order = sortperm(usage; rev=true, alg=Base.Sort.MergeSort)
    theta, omega = best.theta[:, order], best.omega[:, order]
    membership = [argmax(row) for row in eachrow(omega)]
    sizes = [count(==(r), membership) for r in 1:groups]
    return GoMResult(samples, species, ["region_$r" for r in 1:groups], omega, theta,
        groups, membership, sizes, best.likelihood, best.objective, alpha,
        best.iterations, best.converged, best.reason, best.residual, best.history,
        starts, limit, tolerance, accelerate, best_start, objectives, iterations, converged, reasons, residuals)
end

function fitgom(data::CommunityData; kwargs...)
    _validate_community_data(data)
    _community_aggregation_keywords(kwargs)
    return fitgom(data.matrix; sample_labels=data.sample_labels,
        species_labels=data.species_labels, kwargs...)
end

fitgom(community::AbstractMatrix; kwargs...) =
    throw(ArgumentError("community matrix must contain real numbers without missing values"))

"""Pass a fixed-k GoM fit's species probabilities and labels to KL indicators."""
function indicators(fit::GoMResult; kwargs...)
    return indicators(fit.theta; species_labels=fit.species_labels,
        region_labels=fit.region_labels, kwargs...)
end
