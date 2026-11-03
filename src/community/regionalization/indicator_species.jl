"""Ranked species per region, with missing entries when too few are eligible."""
struct IndicatorSpeciesResult
    species_labels::Vector{String}
    region_labels::Vector{String}
    indices::Matrix{Union{Missing, Int}}
    scores::Matrix{Union{Missing, Float64}}
    top_features::Matrix{Union{Missing, String}}
    counts::Vector{Int}
    method::Symbol
    comparison::Symbol
    shared::Bool
end

@inline function _indicator_poisson_kl(p::Float64, q::Float64)
    p == q && return 0.0
    p == 0 && return q
    q == 0 && return Inf
    relative = (p - q) / q
    if abs(relative) < 1e-4
        # (1+r)log(1+r)-r, with no cancellation near equality.
        factor = 0.5 + relative * (-1/6 + relative * (1/12 + relative * (-1/20 + relative/30)))
        return (q * relative) * relative * factor
    end
    logratio = abs(relative) < 0.5 ? log1p(relative) : log(p) - log(q)
    return max(0.0, p * logratio + (q - p))
end

"""
    indicators(theta; species_labels=nothing, region_labels=nothing,
               top_indicators=5, method=:poisson, comparison=:min, shared=false)

Rank species (rows) for each region (columns) by directed KL divergence from
that region to the other regions. Poisson inputs are finite nonnegative
intensities; Bernoulli inputs are probabilities in [0,1]. Columns are not
renormalized. `comparison` aggregates across other regions by `:min` or `:max`.

With `shared=false`, a species is eligible only for its first maximum column.
Scores tie in input species order. Zeros and ones use the continuous KL limits,
including infinite divergence against zero-probability events. Missing ranks
are padded with `missing`; arrays have regions × requested-ranks shape.
Uses O(species × regions) input storage plus one species score vector, avoiding
a species × regions × regions tensor. At least two regions are required.
"""
function indicators(theta::AbstractMatrix{<:Real}; species_labels=nothing, region_labels=nothing,
    top_indicators=5, method::Symbol=:poisson, comparison::Symbol=:min, shared::Bool=false,
)
    n, regions = size(theta)
    n > 0 && regions >= 2 || throw(ArgumentError("theta needs species rows and at least two regions"))
    method in (:poisson, :bernoulli) || throw(ArgumentError("method must be :poisson or :bernoulli"))
    comparison in (:min, :max) || throw(ArgumentError("comparison must be :min or :max"))
    top_indicators isa Integer && !(top_indicators isa Bool) && 1 <= top_indicators <= typemax(Int) ||
        throw(ArgumentError("top_indicators must be a positive Int-sized integer"))
    species = _validated_labels(species_labels, n, "species")
    labels = _validated_labels(region_labels, regions, "region")
    data = Matrix{Float64}(undef, n, regions)
    for j in 1:regions, i in 1:n
        value = theta[i, j]
        isfinite(value) && value >= 0 && (method !== :bernoulli || value <= 1) ||
            throw(ArgumentError("theta values are outside the chosen KL model's domain"))
        converted = Float64(value)
        isfinite(converted) && (value == 0 || converted > 0) &&
            (method !== :bernoulli || value == 1 || converted < 1) ||
            throw(ArgumentError("theta value cannot be represented without losing its probability boundary"))
        data[i, j] = converted
    end
    dominant = [argmax(@view data[i, :]) for i in 1:n]
    indices = Matrix{Union{Missing, Int}}(missing, regions, Int(top_indicators))
    scores = Matrix{Union{Missing, Float64}}(missing, regions, Int(top_indicators))
    features = Matrix{Union{Missing, String}}(missing, regions, Int(top_indicators))
    counts = zeros(Int, regions)
    work = Vector{Float64}(undef, n)
    for region in 1:regions
        for i in 1:n
            score = comparison === :min ? Inf : 0.0
            p = data[i, region]
            for other in 1:regions
                region == other && continue
                q = data[i, other]
                divergence = _indicator_poisson_kl(p, q)
                method === :bernoulli && (divergence += _indicator_poisson_kl(1 - p, 1 - q))
                score = comparison === :min ? min(score, divergence) : max(score, divergence)
            end
            work[i] = score
        end
        eligible = shared ? collect(1:n) : findall(==(region), dominant)
        sort!(eligible; by=i -> (-work[i], i))
        counts[region] = min(length(eligible), Int(top_indicators))
        for rank in 1:counts[region]
            i = eligible[rank]
            indices[region, rank] = i
            scores[region, rank] = work[i]
            features[region, rank] = species[i]
        end
    end
    return IndicatorSpeciesResult(species, labels, indices, scores, features, counts,
        method, comparison, shared)
end
