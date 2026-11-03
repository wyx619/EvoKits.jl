"""Value classes in input order, exact upper bounds, and optional Jenks objective."""
struct ValueClassificationResult
    sample_labels::Vector{String}
    values::Vector{Float64}
    classes::Vector{Union{Missing, Int}}
    breaks::Vector{Float64}
    upper_bounds::Vector{Float64}
    requested_k::Int
    actual_k::Int
    method::Symbol
    objective::Union{Nothing, Float64}
    normalized_objective::Union{Nothing, Float64}
    objective_scale::Float64
    reason::Symbol
end

# Exact weighted 1-D partitioning; equal values are never split across classes.
function _jenks_partition(sorted, k, max_work)
    values, weights = Float64[], Int[]
    for value in sorted
        if isempty(values) || value != last(values)
            push!(values, value)
            push!(weights, 1)
        else
            weights[end] += 1
        end
    end
    m = length(values)
    groups = min(k, m)
    BigInt(groups) * m * (m + 1) ÷ 2 <= max_work ||
        throw(ArgumentError("Jenks exact dynamic program exceeds max_work; reduce distinct values or raise the budget"))
    span = last(values) - first(values)
    scale = span > 0 && isfinite(span) ? span : max(abs(first(values)), abs(last(values)), 1.0)
    normalized = isfinite(span) ? (values .- first(values)) ./ scale :
        values ./ scale .- first(values) / scale
    costs = fill(Inf, groups, m)
    starts = zeros(Int, groups, m)
    for last_index in 1:m
        mean, sse, count = 0.0, 0.0, 0
        for first_index in last_index:-1:1
            weight = weights[first_index]
            next_count = count + weight
            delta = normalized[first_index] - mean
            next_mean = mean + delta * (weight / next_count)
            sse += weight * delta * (normalized[first_index] - next_mean)
            mean, count = next_mean, next_count
            if first_index == 1
                costs[1, last_index], starts[1, last_index] = max(sse, 0.0), 1
            else
                for group in 2:min(groups, first_index)
                    candidate = costs[group - 1, first_index - 1] + max(sse, 0.0)
                    if candidate < costs[group, last_index]
                        costs[group, last_index], starts[group, last_index] = candidate, first_index
                    end
                end
            end
        end
    end
    upper = zeros(Float64, groups)
    last_index = m
    for group in groups:-1:1
        upper[group] = values[last_index]
        last_index = starts[group, last_index] - 1
    end
    normalized_sse = costs[groups, m]
    # May overflow in original squared units; the normalized objective stays finite.
    root_sse = sqrt(normalized_sse) * scale
    return upper, root_sse * root_sse, normalized_sse, scale
end

"""
    classify_values(values; k=5, method=:jenks, sample_labels=nothing,
                    invalid=:error, max_work=50_000_000)

Classify a metric by exact Jenks natural breaks, equal widths or type-7 quantiles.
Jenks minimizes total within-class squared deviations using sorted unique values
with multiplicity weights and Welford interval variances. Its time is O(k m²)
and workspace O(k m), with m distinct values; `max_work` bounds that workload.
Translation/scaling avoids raw squared-sum cancellation. Original-unit SSE can
overflow; `normalized_objective` and `objective_scale` retain the scaled result.

Zero, negative and fractional values are retained. Tied values share a class.
Jenks caps classes at the distinct-value count; equal/quantile methods cap bins
at the observation count and deduplicate upper bounds. `actual_k` counts those
bounds; equal/quantile bins can be empty. A reduced count gives `:reduced_classes`. Invalid
missing/NaN positions remain missing under `invalid=:skip`; Inf always fails.
Classes include their upper boundary. `breaks` includes the minimum then all
upper bounds (the first two may coincide for a singleton minimum class).
Equal-cost Jenks splits prefer the latest final split, recursively. At least
one valid observation is required. Inputs are not modified.
"""
function classify_values(values::AbstractVector; k=5, method::Symbol=:jenks,
    sample_labels=nothing, invalid::Symbol=:error, max_work=50_000_000,
)
    k isa Integer && !(k isa Bool) && 1 <= k <= typemax(Int) || throw(ArgumentError("k must be a positive Int-sized integer"))
    max_work isa Integer && !(max_work isa Bool) && max_work > 0 || throw(ArgumentError("max_work must be a positive integer"))
    method in (:jenks, :equal, :quantile) || throw(ArgumentError("method must be :jenks, :equal or :quantile"))
    data, valid = _metric_values(values; invalid)
    labels = _community_data_labels(sample_labels, length(data), "sample")
    sorted = sort(data[valid])
    isempty(sorted) && throw(ArgumentError("classification needs at least one valid value"))
    objective, normalized, scale = nothing, nothing, 1.0
    requested = Int(k)
    if method === :jenks
        upper, objective, normalized, scale = _jenks_partition(sorted, requested, max_work)
    else
        # Avoid creating more empty bins than observations.
        groups = min(requested, length(sorted))
        upper = method === :equal ?
            [_metric_lerp(first(sorted), last(sorted), i / groups) for i in 1:groups] :
            [_metric_quantile(sorted, i / groups) for i in 1:groups]
        unique!(upper)
    end
    classes = Vector{Union{Missing, Int}}(missing, length(data))
    for i in eachindex(data)
        valid[i] && (classes[i] = searchsortedfirst(upper, data[i]))
    end
    return ValueClassificationResult(labels, data, classes, vcat(first(sorted), upper), upper,
        requested, length(upper), method, objective, normalized, scale,
        length(upper) == requested ? :ok : :reduced_classes)
end
