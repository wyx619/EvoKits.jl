"""
    MeanDistResult

Distance average returned by [`mean_dist`](@ref). `distances` and `counts` are
symmetric matrices, or vectors in `pairs` order. `counts` records the number of
finite observations for each pair; `replicates` counts all consumed inputs.
`sample_labels` fixes the output site order. `component` is the selected beta
component, or :distance for matrix inputs; `index_family` is nothing for plain
distances. Results own their arrays.
"""
struct MeanDistResult{N}
    sample_labels::Vector{String}
    distances::Array{Float64, N}
    counts::Array{Int, N}
    pairs::Union{Nothing, Vector{Tuple{Int, Int}}}
    replicates::Int
    component::Symbol
    index_family::Union{Nothing, Symbol}
end

function _mean_dist_labels(labels)
    labels isa AbstractVector || throw(ArgumentError("sample_labels must be a vector"))
    return _community_data_labels(labels, length(labels), "sample")
end

function _mean_dist_input(input::Union{BetaDissResult, PhylobetaResult}, component, labels)
    _validate_community_beta_family(input.index_family)
    return (distances = getproperty(input, component),
        labels = _mean_dist_labels(input.sample_labels), pairs = input.pairs,
        kind = input isa PhylobetaResult ? :phylogenetic : :taxonomic,
        family = input.index_family)
end

function _mean_dist_input(input::AbstractMatrix{<:Real}, component, labels)
    return (distances = input,
        labels = _community_data_labels(labels, size(input, 1), "sample"),
        pairs = nothing, kind = :distance, family = nothing)
end

function _mean_dist_input(input::NamedTuple, component, labels)
    haskey(input, :distances) && haskey(input, :sample_labels) ||
        throw(ArgumentError("labeled inputs require distances and sample_labels"))
    return (distances = input.distances, labels = _mean_dist_labels(input.sample_labels),
        pairs = get(input, :pairs, nothing), kind = :distance, family = nothing)
end

_mean_dist_input(input, component, labels) = throw(ArgumentError(
    "inputs must contain beta results, distance matrices, or labeled distance named tuples"))

function _mean_dist_value(value)
    value isa Real || throw(ArgumentError("distances must be real numbers or NaN"))
    isnan(value) && return NaN
    isfinite(value) && value >= 0 ||
        throw(ArgumentError("distances must be non-negative and finite, or NaN"))
    converted = Float64(value)
    isfinite(converted) && (value == 0 || converted > 0) || throw(ArgumentError(
        "distances must be representable as finite Float64 without underflow to zero"))
    return converted
end

function _mean_dist_matrix_value(distances, i, j)
    value, reverse = distances[i, j], distances[j, i]
    converted = _mean_dist_value(value)
    _mean_dist_value(reverse)
    (value == reverse || (isnan(value) && isnan(reverse))) ||
        throw(ArgumentError("distance matrices must be exactly symmetric, including NaN locations"))
    i != j || isnan(converted) || converted == 0 ||
        throw(ArgumentError("distance diagonals must be zero or NaN"))
    return converted
end

# Scale a compensated running sum by the largest observed distance. The sum
# cannot overflow for an Int-sized replicate count, and tiny original units
# are retained until the final division (including subnormal distances).
@inline function _mean_dist_add!(state, index, value, nan_policy)
    if isnan(value)
        nan_policy === :propagate && (state.sums[index] = NaN)
        return nothing
    end
    state.counts[index] += 1
    old_sum = state.sums[index]
    isnan(old_sum) && return nothing
    scale = state.scales[index]
    if value > scale
        factor = scale / value
        old_sum *= factor
        state.corrections[index] *= factor
        state.scales[index] = scale = value
    end
    term = scale == 0 ? 0.0 : value / scale
    updated = old_sum + term
    state.corrections[index] += old_sum >= term ?
        (old_sum - updated) + term : (term - updated) + old_sum
    state.sums[index] = updated
    return nothing
end

@inline function _mean_dist_index(i, j, n)
    lo, hi = minmax(i, j)
    return hi + (lo - 1) * n
end

function _mean_dist_full!(state, distances, mapping, selected, nan_policy, n)
    m = length(mapping)
    # Validate the complete supplied matrix, also when selecting output pairs.
    for j in 1:m, i in j:m
        value = _mean_dist_matrix_value(distances, i, j)
        selected === nothing && _mean_dist_add!(state,
            _mean_dist_index(mapping[i], mapping[j], n), value, nan_policy)
    end
    if selected !== nothing
        inverse = zeros(Int, n)
        for (i, global_index) in enumerate(mapping)
            inverse[global_index] = i
        end
        for (index, (a, b)) in enumerate(selected)
            i, j = inverse[a], inverse[b]
            i > 0 && j > 0 || continue
            _mean_dist_add!(state, index, _mean_dist_value(distances[i, j]), nan_policy)
        end
    end
    return nothing
end

function _mean_dist_pairs!(state, distances, input_pairs, mapping, selected, nan_policy)
    local_pairs = _community_pairs(input_pairs, length(mapping))
    length(local_pairs) == length(distances) ||
        throw(ArgumentError("distance vectors must match the input pairs length"))
    values = Dict{Tuple{Int, Int}, Float64}()
    for (index, (i, j)) in enumerate(local_pairs)
        value = _mean_dist_value(distances[index])
        i != j || isnan(value) || value == 0 ||
            throw(ArgumentError("self-pair distances must be zero or NaN"))
        key = minmax(mapping[i], mapping[j])
        if haskey(values, key)
            previous = values[key]
            (previous == value || (isnan(previous) && isnan(value))) ||
                throw(ArgumentError("repeated or reversed input pairs must have equal distances"))
        else
            values[key] = value
        end
    end
    # One observation per input and unordered pair: repeated source pairs do
    # not give a posterior draw extra weight. Output repetitions are preserved.
    for (index, (i, j)) in enumerate(selected)
        key = minmax(i, j)
        haskey(values, key) || continue
        _mean_dist_add!(state, index, values[key], nan_policy)
    end
    return nothing
end

function _mean_dist_finish!(state, selected, n)
    indices = selected === nothing ?
        (_mean_dist_index(i, j, n) for j in 1:n for i in j:n) : eachindex(state.sums)
    for index in indices
        count, sum = state.counts[index], state.sums[index]
        state.sums[index] = count == 0 || isnan(sum) ? NaN :
            state.scales[index] * clamp((sum + state.corrections[index]) / count, 0.0, 1.0)
    end
    if selected === nothing
        for j in 1:n, i in (j + 1):n
            state.sums[j, i] = state.sums[i, j]
            state.counts[j, i] = state.counts[i, j]
        end
    end
    return nothing
end

"""
    mean_dist(inputs; component=:total, sample_labels=nothing, pairs=nothing,
              nan_policy=:skip)

Average distance observations across posterior trees or other replicates.
`inputs` is a nonempty iterable consumed once; a generator may compute one
`PhylobetaResult` or `BetaDissResult` at a time. Select :total, :turnover or
:nestedness. All beta inputs must have the same kind and index family.

Alternatively supply matrices, or `(distances=D, sample_labels=labels)` named
tuples; a labeled vector also needs `pairs`. Plain distances cannot be mixed
with beta results. Bare matrices use the common output label order and must
contain its complete domain. Matrix results record component=:distance.

Output labels preserve the first input's order, or the explicit `sample_labels`
order. Subsequent labeled inputs may be reordered or contain subsets. Labels
match exactly; unknown, duplicate or empty labels fail. To include sites not
in the first input, supply the complete output domain explicitly.

The default :skip ignores each NaN independently and counts only finite values.
Pairs with no finite observations return NaN with count zero, including empty
site diagonals. :propagate returns NaN if any observed value is NaN, while
`counts` still reports finite observations. Absent sites/pairs are unobserved,
not NaN observations. Infinities and negative values fail. Finite diagonals
must be zero; full inputs must be exactly symmetric, including NaN locations.
Each beta component has its own finite-value mask, so averages of separately
selected components need not sum to the average total.

Supply output `pairs` as indices into the output labels to return vectors,
preserving repeated, reversed and self pairs. Pair-vector inputs require this
explicit option; their own pair indices refer to their own labels. Repeated
source pairs must agree and contribute once per replicate. Selected output
uses O(sites + requested pairs + current input pairs) auxiliary space, with no
full matrix constructed. Full output and its work buffers use O(sites^2).
Complete supplied matrices are validated even for selected output. Inputs are
not retained or mutated; a scaled compensated sum avoids overflow in averaging.

Pass a complete, finite result with a zero diagonal directly to `phyloregion`,
`prepare_phyloregion`, `select_linkage`, or `optimal_phyloregion`. Their existing
rules reject undefined sites and selected-pair results; no site is dropped or
distance imputed. This summarizes distances, not partition uncertainty or
posterior support for individual regions.
"""
function mean_dist(inputs; component::Symbol = :total, sample_labels = nothing,
    pairs = nothing, nan_policy::Symbol = :skip,
)
    component in (:total, :turnover, :nestedness) ||
        throw(ArgumentError("component must be :total, :turnover or :nestedness"))
    nan_policy in (:skip, :propagate) ||
        throw(ArgumentError("nan_policy must be :skip or :propagate"))
    labels = sample_labels === nothing ? nothing : _mean_dist_labels(sample_labels)
    selected, state, lookup = nothing, nothing, nothing
    kind, family, replicates = nothing, nothing, 0
    for input in inputs
        source = _mean_dist_input(input, component, labels)
        # Reject invalid layouts before allocating a possibly quadratic output.
        if source.pairs === nothing
            source.distances isa AbstractMatrix{<:Real} ||
                throw(ArgumentError("inputs without pairs require real distance matrices"))
            m = length(source.labels)
            size(source.distances) == (m, m) || throw(ArgumentError(
                "distance matrices must be square and match sample_labels"))
        else
            pairs !== nothing || throw(ArgumentError(
                "pair-vector inputs require explicit output pairs; no full matrix is inferred"))
            source.distances isa AbstractVector{<:Real} ||
                throw(ArgumentError("inputs with pairs require real distance vectors"))
        end
        Base.require_one_based_indexing(source.distances)
        if replicates == 0
            labels = labels === nothing ? copy(source.labels) : labels
            selected = _community_pairs(pairs, length(labels))
            shape = selected === nothing ? (length(labels), length(labels)) : (length(selected),)
            state = (sums = zeros(Float64, shape), counts = zeros(Int, shape),
                scales = zeros(Float64, shape), corrections = zeros(Float64, shape))
            lookup = Dict(label => i for (i, label) in enumerate(labels))
            kind, family = source.kind, source.family
        else
            (source.kind, source.family) == (kind, family) || throw(ArgumentError(
                "inputs must have the same distance kind and beta index family"))
        end
        mapping = [get(lookup, label, 0) for label in source.labels]
        all(>(0), mapping) || throw(ArgumentError(
            "input sample labels must belong to the output label domain"))
        if source.pairs === nothing
            _mean_dist_full!(state, source.distances, mapping, selected, nan_policy, length(labels))
        else
            _mean_dist_pairs!(state, source.distances, source.pairs, mapping, selected, nan_policy)
        end
        replicates += 1
    end
    replicates > 0 || throw(ArgumentError("inputs must not be empty"))
    _mean_dist_finish!(state, selected, length(labels))
    return MeanDistResult(labels, state.sums, state.counts, selected, replicates,
        kind === :distance ? :distance : component, family)
end

function _mean_dist_region_matrix(result::MeanDistResult)
    result.pairs === nothing && result.distances isa AbstractMatrix ||
        throw(ArgumentError("regionalization requires a full mean distance matrix, not selected pairs"))
    return result.distances
end

function phyloregion(result::MeanDistResult; k, method::Symbol = :average)
    return _phyloregion(_mean_dist_region_matrix(result), k, result.sample_labels,
        method, result.component)
end

function prepare_phyloregion(result::MeanDistResult; method::Symbol = :average)
    return _prepare_phyloregion(_mean_dist_region_matrix(result), result.sample_labels,
        method, result.component)
end

function select_linkage(result::MeanDistResult; methods = _REGION_METHODS)
    distances = _mean_dist_region_matrix(result)
    return _select_linkage(_region_distance_buffer(distances), size(distances, 1),
        result.component, methods)
end

function optimal_phyloregion(result::MeanDistResult; method::Symbol = :average, kwargs...)
    return optimal_phyloregion(prepare_phyloregion(result; method); kwargs...)
end
