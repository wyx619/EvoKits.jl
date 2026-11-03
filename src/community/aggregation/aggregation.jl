"""Per-site numeric trait summaries returned by [`map_trait`](@ref)."""
struct CommunityTraitResult
    sample_labels::Vector{String}
    values::Vector{Float64}
    richness::Vector{Int}
    matched_richness::Vector{Int}
    weight_sum::Vector{Float64}
    statistic::Symbol
    abundance_weighted::Bool
    corrected::Bool
end

# All kernels visit columns in input order, including CSC columns with explicit
# zeros. No dense community, long occurrence table or membership matrix is built.
function _community_each_positive(f, matrix::AbstractMatrix{<:Real})
    for col in axes(matrix, 2), row in axes(matrix, 1)
        value = matrix[row, col]
        value > 0 && f(row, col, value)
    end
    return nothing
end

function _community_each_positive(f, matrix::SparseMatrixCSC{<:Real})
    for col in axes(matrix, 2), entry in nzrange(matrix, col)
        value = nonzeros(matrix)[entry]
        value > 0 && f(rowvals(matrix)[entry], col, value)
    end
    return nothing
end

function _community_mapping(identifiers, values)
    length(identifiers) == length(values) || throw(ArgumentError("species and values must have equal lengths"))
    labels = _validated_labels(identifiers, length(identifiers), "species")
    return Dict(label => value for (label, value) in zip(labels, values))
end

function _community_mapping(mapping::AbstractDict)
    labels = collect(keys(mapping))
    return _community_mapping(labels, [mapping[label] for label in labels])
end

function _community_mapping_values(mapping::AbstractVector, labels)
    Base.require_one_based_indexing(mapping)
    length(mapping) == length(labels) || throw(ArgumentError("mapping vector must have one value per community species"))
    return mapping
end

function _community_mapping_values(mapping::AbstractDict, labels)
    normalized = _community_mapping(mapping)
    return [get(normalized, label, missing) for label in labels]
end

function _community_mapping_values(mapping, labels)
    throw(ArgumentError("mapping must be a vector, a species-keyed dictionary, or a DataFrame"))
end

function _community_aggregation_input(matrix::AbstractMatrix{<:Real}, species_labels, sample_labels)
    _validate_community_values(matrix)
    species = _validated_labels(species_labels, size(matrix, 2), "species")
    samples = _validated_labels(sample_labels, size(matrix, 1), "sample")
    return species, samples
end

function _community_missing_mapping(policy, name)
    policy in (:error, :skip) || throw(ArgumentError("$name must be :error or :skip"))
end

function _community_aggregation_keywords(kwargs)
    _community_data_keywords(kwargs)
    haskey(kwargs, :species_labels) && throw(ArgumentError("CommunityData owns its labels; set species_labels when constructing the data"))
end

function _community_trait_values(mapping, species)
    raw = _community_mapping_values(mapping, species)
    values = zeros(Float64, length(raw))
    valid = falses(length(raw))
    for (col, value) in enumerate(raw)
        ismissing(value) && continue
        value isa Real && isfinite(value) || throw(ArgumentError("traits must be finite real numbers or missing"))
        converted = Float64(value)
        isfinite(converted) || throw(ArgumentError("traits must be representable as finite Float64 values"))
        values[col], valid[col] = converted, true
    end
    return values, valid
end

"""
    map_trait(community, traits; statistic=:sum, abundance_weighted=false,
              corrected=true, missing_traits=:error, species_labels=nothing,
              sample_labels=nothing)
    map_trait(data::CommunityData, traits; kwargs...)
    map_trait(community, table::AbstractDataFrame; species=:species, trait=:trait, kwargs...)

Summarize one numeric species trait per site (`:sum`, `:mean`, or `:sd`). Traits
are a vector in community column order, a dictionary keyed by species labels,
or a table selected by `species` and `trait`. Names match exactly; duplicate
table species are errors. Extra mapping species are ignored. A missing mapping
or `missing` value is an error for a present species, or is omitted with
`missing_traits=:skip`. Unobserved species need not have a trait. Nonfinite and
nonnumeric matched trait values are always errors.

Positive community values denote presence. With `abundance_weighted=true`, use
their abundances as weights: sum is sum(w*x), mean is sum(w*x)/sum(w), and SD
uses weighted squared deviations. `corrected=true` uses denominator
`sum(w) - sum(w^2)/sum(w)` (sample SD for unit weights); `corrected=false` uses
sum(w). Weights are reliability weights, so weighted mean/SD are invariant to
common abundance scaling. Corrected SD requires at least two matched species.

Returns `CommunityTraitResult` with `values`, input `richness`, used
`matched_richness`, used `weight_sum`, and site labels. Empty/unmatched sites
have zero sum, NaN mean/SD, and zero matched richness/weight. Calculations use
Float64 and reject overflowing weights or summaries. Dense and CSC matrices
are supported; CSC is never expanded. Inputs and result arrays are independent.
"""
function map_trait(community::AbstractMatrix{<:Real}, traits;
    statistic::Symbol = :sum, abundance_weighted::Bool = false,
    corrected::Bool = true, missing_traits::Symbol = :error,
    species_labels = nothing, sample_labels = nothing,
)
    statistic in (:sum, :mean, :sd) || throw(ArgumentError("statistic must be :sum, :mean or :sd"))
    _community_missing_mapping(missing_traits, "missing_traits")
    species, samples = _community_aggregation_input(community, species_labels, sample_labels)
    trait_values, valid = _community_trait_values(traits, species)
    n = size(community, 1)
    richness, matched = zeros(Int, n), zeros(Int, n)
    weights, values, means, m2, correction = (zeros(n) for _ in 1:5)
    anchors = statistic === :sd ? zeros(n) : Float64[]
    _community_each_positive(community) do row, col, abundance
        richness[row] += 1
        if !valid[col]
            missing_traits === :error && throw(ArgumentError("missing trait for present species: $(species[col])"))
            return
        end
        weight = abundance_weighted ? Float64(abundance) : 1.0
        weight > 0 && isfinite(weight) || throw(ArgumentError("abundances must be representable as positive finite Float64 weights"))
        previous = weights[row]
        total = previous + weight
        isfinite(total) || throw(ArgumentError("trait weight sum overflowed"))
        x = trait_values[col]
        if statistic === :sum
            term = weight * x
            total_value = values[row] + term
            # Neumaier summation retains small signed traits during cancellation.
            means[row] += abs(values[row]) >= abs(term) ?
                (values[row] - total_value) + term : (term - total_value) + values[row]
            values[row] = total_value
        elseif previous == 0
            means[row] = statistic === :sd ? 0.0 : x
            statistic === :sd && (anchors[row] = x)
        else
            # Center SD at an actual trait so large shared offsets do not
            # erase a small between-species spread in the running mean.
            statistic === :sd && (x -= anchors[row])
            delta = x - means[row]
            ratio = weight / total
            means[row] = isfinite(delta) ? means[row] + delta * ratio :
                means[row] * (previous / total) + x * ratio
            if statistic === :sd
                m2[row] += (previous * ratio) * delta * delta
                # Accumulate W - sum(w^2)/W directly, avoiding cancellation
                # when one species has almost all of a site's weight.
                correction[row] = correction[row] * (previous / total) + 2 * (weight * (previous / total))
            end
        end
        weights[row] = total
        matched[row] += 1
    end
    for row in 1:n
        if statistic === :sum
            values[row] += means[row]
        elseif statistic === :mean
            values[row] = matched[row] > 0 ? means[row] : NaN
        elseif statistic === :sd
            denominator = corrected ? correction[row] : weights[row]
            values[row] = matched[row] > 0 && denominator > 0 ? sqrt(max(0.0, m2[row] / denominator)) : NaN
        end
        (matched[row] == 0 || (statistic === :sd && corrected && matched[row] == 1)) && continue
        isfinite(values[row]) || throw(ArgumentError("trait summary overflowed; rescale traits or abundances"))
    end
    return CommunityTraitResult(samples, values, richness, matched, weights,
        statistic, abundance_weighted, corrected)
end

function map_trait(community::AbstractMatrix{<:Real}, table::AbstractDataFrame;
    species = :species, trait = :trait, kwargs...,
)
    mapping = _community_mapping(_community_table_column(table, species), _community_table_column(table, trait))
    return map_trait(community, mapping; kwargs...)
end

function map_trait(data::CommunityData, traits; kwargs...)
    _community_aggregation_keywords(kwargs)
    _validate_community_data(data)
    return map_trait(data.matrix, traits; species_labels = data.species_labels,
        sample_labels = data.sample_labels, kwargs...)
end

function _community_group_indices(mapping, species, group_labels)
    raw = _community_mapping_values(mapping, species)
    groups = [ismissing(value) ? missing : _label_identifier(value, "group identifiers") for value in raw]
    domain = group_labels === nothing ? String[value for value in unique(groups) if !ismissing(value)] :
        _validated_labels(group_labels, length(group_labels), "group")
    positions = Dict(label => i for (i, label) in enumerate(domain))
    indices = zeros(Int, length(groups))
    for (col, group) in enumerate(groups)
        ismissing(group) && continue
        indices[col] = get(positions, group, 0)
        indices[col] > 0 || throw(ArgumentError("mapped group $group is absent from group_labels"))
    end
    return indices, domain
end

function _community_normalize_groups(matrix::SparseMatrixCSC{T}) where {T}
    maximums = zeros(T, size(matrix, 1))
    for entry in eachindex(nonzeros(matrix))
        row = rowvals(matrix)[entry]
        maximums[row] = max(maximums[row], nonzeros(matrix)[entry])
    end
    values = [Float64(value / maximums[row]) for (row, value) in zip(rowvals(matrix), nonzeros(matrix))]
    totals = zeros(Float64, size(matrix, 1))
    for entry in eachindex(values)
        totals[rowvals(matrix)[entry]] += values[entry]
    end
    for entry in eachindex(values)
        values[entry] /= totals[rowvals(matrix)[entry]]
    end
    return SparseMatrixCSC(size(matrix, 1), size(matrix, 2), copy(matrix.colptr), copy(rowvals(matrix)), values)
end

"""
    aggregate_species(community, groups; abundance_weighted=false, normalize=false,
                      missing_groups=:error, group_labels=nothing,
                      species_labels=nothing, sample_labels=nothing)
    aggregate_species(data::CommunityData, groups; kwargs...)
    aggregate_species(community, table::AbstractDataFrame; species=:species, group=:group, kwargs...)

Aggregate species columns into predefined groups. `groups` is a vector in
community column order, a species-keyed dictionary, or a two-column table.
Group identifiers are nonempty strings or integers; each species maps to at
most one group. Table species must be unique. Missing mappings/`missing` values
are errors for present species (`missing_groups=:error`) or can be omitted
(`:skip`). Unobserved species may be unmapped. Extra mapping species are ignored.

By default count each present species once per site. With
`abundance_weighted=true`, sum original abundance instead, preserving numeric
type (Bool becomes Int); overflowing sums fail. `normalize=true` returns
Float64 proportions of retained counts/abundances per site. Empty rows remain
zero, including with normalization.

Returns a new `CommunityData` with a CSC matrix, original site order/labels, and
group identifiers as `species_labels`. Groups default to first occurrence in
community species-column order, including mapped zero columns. Supply complete
`group_labels` to set order and retain extra empty groups. Dense and CSC inputs
are supported without building a dense community or group-membership matrix.
No input is mutated. Clustering and tree cutting are separate from this operation.
"""
function aggregate_species(community::AbstractMatrix{<:Real}, groups;
    abundance_weighted::Bool = false, normalize::Bool = false,
    missing_groups::Symbol = :error, group_labels = nothing,
    species_labels = nothing, sample_labels = nothing,
)
    _community_missing_mapping(missing_groups, "missing_groups")
    species, samples = _community_aggregation_input(community, species_labels, sample_labels)
    indices, labels = _community_group_indices(groups, species, group_labels)
    T = abundance_weighted && eltype(community) !== Bool ? eltype(community) : Int
    rows, cols, values = Int[], Int[], T[]
    _community_each_positive(community) do row, col, abundance
        index = indices[col]
        if index == 0
            missing_groups === :error && throw(ArgumentError("missing group for present species: $(species[col])"))
            return
        end
        push!(rows, row)
        push!(cols, index)
        push!(values, abundance_weighted ? abundance : one(T))
    end
    matrix = sparse(rows, cols, values, size(community, 1), length(labels), _community_duplicate_sum)
    dropzeros!(matrix)
    result = normalize ? _community_normalize_groups(matrix) : matrix
    return CommunityData(result, samples, labels, Val(:validated))
end

function aggregate_species(community::AbstractMatrix{<:Real}, table::AbstractDataFrame;
    species = :species, group = :group, kwargs...,
)
    mapping = _community_mapping(_community_table_column(table, species), _community_table_column(table, group))
    return aggregate_species(community, mapping; kwargs...)
end

function aggregate_species(data::CommunityData, groups; kwargs...)
    _community_aggregation_keywords(kwargs)
    _validate_community_data(data)
    return aggregate_species(data.matrix, groups; species_labels = data.species_labels,
        sample_labels = data.sample_labels, kwargs...)
end
