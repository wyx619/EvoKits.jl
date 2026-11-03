# Shared validation for vector-based conservation summaries and classifications.
function _metric_values(values::AbstractVector; invalid::Symbol=:error)
    invalid in (:error, :skip) || throw(ArgumentError("invalid must be :error or :skip"))
    result = fill(NaN, length(values))
    valid = falses(length(values))
    for (i, value) in enumerate(values)
        if ismissing(value) || (value isa Real && isnan(value))
            invalid === :skip || throw(ArgumentError("metric values contain missing or NaN"))
            continue
        end
        value isa Real && isfinite(value) || throw(ArgumentError("metric values must be finite real numbers"))
        converted = Float64(value)
        isfinite(converted) && (value == 0 || converted != 0) ||
            throw(ArgumentError("metric value is not representable as finite Float64"))
        result[i], valid[i] = converted, true
    end
    return result, valid
end

@inline function _metric_lerp(lo, hi, fraction)
    fraction == 0 && return lo
    fraction == 1 && return hi
    span = hi - lo
    return isfinite(span) ? lo + span * fraction : lo * (1 - fraction) + hi * fraction
end

# R type-7 quantiles, with safe interpolation across very large opposite signs.
function _metric_quantile(sorted, probability)
    rank = (length(sorted) - 1) * probability
    lower = floor(Int, rank) + 1
    lower == length(sorted) && return last(sorted)
    return _metric_lerp(sorted[lower], sorted[lower + 1], rank - floor(rank))
end

"""A labeled tail selection; invalid records keep a missing selection flag."""
struct MetricSelectionResult
    sample_labels::Vector{String}
    values::Vector{Float64}
    selected::Vector{Union{Missing, Bool}}
    threshold::Float64
    prob::Float64
    tail::Symbol
    valid_count::Int
    reference_count::Int
end

function _metric_tail(values; sample_labels=nothing, reference=nothing,
    prob=2.5, invalid::Symbol=:error, tail,
)
    prob isa Real && !(prob isa Bool) && isfinite(prob) && 0 < prob <= 100 &&
        0 < Float64(prob) <= 100 || throw(ArgumentError("prob is a percentage in (0,100]"))
    data, valid = _metric_values(values; invalid)
    labels = _community_data_labels(sample_labels, length(data), "sample")
    refs, reference_valid = reference === nothing ? (data, valid) : _metric_values(reference; invalid)
    sorted = sort(refs[reference_valid])
    isempty(sorted) && throw(ArgumentError("the threshold reference has no valid values"))
    probability = Float64(prob) / 100
    threshold = _metric_quantile(sorted, tail === :high ? 1 - probability : probability)
    flags = Vector{Union{Missing, Bool}}(missing, length(data))
    for i in eachindex(data)
        valid[i] && (flags[i] = tail === :high ? data[i] >= threshold : data[i] <= threshold)
    end
    return MetricSelectionResult(labels, data, flags, threshold, Float64(prob), tail,
        count(valid), length(sorted))
end

"""
    hotspots(values; sample_labels=nothing, reference=nothing, prob=2.5, invalid=:error)
    coldspots(values; same keywords...)

Select values at or above the upper, or at or below the lower, type-7 quantile.
`prob` is a tail percentage in (0,100]. Ties are all retained, so the selected
fraction can exceed `prob`. A separate reference vector supplies a common
threshold across datasets. `invalid=:skip` preserves missing/NaN positions
with missing selection flags; Inf is always rejected. A reference must have
at least one valid value. PD, WE and PE results supply values and site labels.
"""
hotspots(values::AbstractVector; kwargs...) = _metric_tail(values; tail=:high, kwargs...)
coldspots(values::AbstractVector; kwargs...) = _metric_tail(values; tail=:low, kwargs...)

for (T, field) in ((PDResult, :pd), (WeightedEndemismResult, :we), (PhyloEndemismResult, :pe))
    @eval begin
        hotspots(result::$T; kwargs...) = hotspots(getfield(result, $(QuoteNode(field)));
            sample_labels=result.sample_labels, kwargs...)
        coldspots(result::$T; kwargs...) = coldspots(getfield(result, $(QuoteNode(field)));
            sample_labels=result.sample_labels, kwargs...)
    end
end

"""Fair-proportion ED and explicitly defined GE scores in complete tree-tip order."""
struct EDGEResult
    species_labels::Vector{String}
    ed::Vector{Float64}
    ge::Vector{Union{Missing, Float64}}
    edge::Vector{Union{Missing, Float64}}
    status::Vector{Symbol}
    ge_mapping::Dict{String, Union{Missing, Float64}}
end

function _edge_mapping(mapping)
    mapping === nothing && return Dict{String, Union{Missing, Float64}}(
        "LC"=>0.001, "NT"=>0.01, "VU"=>0.1, "EN"=>0.67, "CR"=>0.999,
        "EW"=>1.0, "EX"=>1.0, "DD"=>missing)
    mapping isa AbstractDict || throw(ArgumentError("ge_mapping must be a category dictionary"))
    result = Dict{String, Union{Missing, Float64}}()
    for (key, value) in mapping
        key isa Union{Symbol, AbstractString} || throw(ArgumentError("category names must be text"))
        label = String(key)
        !isempty(strip(label)) && !haskey(result, label) || throw(ArgumentError("empty or duplicate category"))
        if ismissing(value)
            result[label] = missing
        else
            value isa Real && !(value isa Bool) && isfinite(value) && value >= 0 ||
                throw(ArgumentError("GE scores must be finite and nonnegative or missing"))
            converted = Float64(value)
            isfinite(converted) && (value == 0 || converted > 0) || throw(ArgumentError("GE score is not representable"))
            result[label] = converted
        end
    end
    return result
end

"""
    edge(tree, categories; tip_labels=tree.tip_labels, ge_mapping=nothing,
         missing_category=:skip)
    edge(tree, table; species_col=:species, category_col=:category, kwargs...)

Compute `log1p(fair_proportion_ED) + GE * log(2)` on the complete input tree.
The default GE mapping matches phyloregion: LC=.001, NT=.01, VU=.1, EN=.67,
CR=.999, EW/EX=1, DD=missing. These are probability-based scores, distinct from
the integer category codes used in some EDGE conventions. Supply a complete
custom `ge_mapping` to choose another convention. This is separate from EDGE2.

Labels match exactly; subsets are accepted and output retains every tree tip.
Missing categories, DD and absent records give missing scores and status under
`:skip`; `:error` rejects them. Unknown species, unknown categories and duplicate
records always fail. Inputs are not modified; branches must be finite/nonnegative.
"""
function edge(tree::CompactTree, categories::AbstractVector; tip_labels=tree.tip_labels,
    ge_mapping=nothing, missing_category::Symbol=:skip,
)
    missing_category in (:skip, :error) || throw(ArgumentError("missing_category must be :skip or :error"))
    _community_data_labels(tree.tip_labels, tree.ntips, "species")
    all(x -> isfinite(x) && x >= 0, tree.edge_length) || throw(ArgumentError("invalid tree branch lengths"))
    labels = _community_data_labels(tip_labels, length(categories), "species")
    positions = Dict(label => i for (i, label) in enumerate(tree.tip_labels))
    mapping = _edge_mapping(ge_mapping)
    ge = Vector{Union{Missing, Float64}}(missing, tree.ntips)
    scores = similar(ge)
    fill!(scores, missing)
    status = fill(:absent_record, tree.ntips)
    for (i, label) in enumerate(labels)
        haskey(positions, label) || throw(ArgumentError("species $label is not in the tree"))
        slot = positions[label]
        category = categories[i]
        if ismissing(category)
            status[slot] = :missing_category
        else
            category isa Union{AbstractString, Symbol} || throw(ArgumentError("categories must be strings, symbols or missing"))
            name = String(category)
            haskey(mapping, name) || throw(ArgumentError("unknown category: $name"))
            ge[slot] = mapping[name]
            status[slot] = ismissing(ge[slot]) ? :unassessed : :ok
        end
    end
    missing_category === :error && any(!=(:ok), status) && throw(ArgumentError("missing or unassessed species categories"))
    ed = evol_distinct(tree; type=:fair_proportion).scores
    all(x -> isfinite(x) && x >= 0, ed) || throw(ArgumentError("ED calculation overflowed"))
    for i in eachindex(ed)
        if status[i] === :ok
            value = log1p(ed[i]) + ge[i] * log(2.0)
            isfinite(value) || throw(ArgumentError("EDGE score overflowed"))
            scores[i] = value
        end
    end
    return EDGEResult(copy(tree.tip_labels), ed, ge, scores, status, mapping)
end

function edge(tree::CompactTree, table::AbstractDataFrame; species_col=:species,
    category_col=:category, kwargs...,
)
    String(species_col) in names(table) && String(category_col) in names(table) ||
        throw(ArgumentError("category table is missing selected columns"))
    return edge(tree, table[!, category_col]; tip_labels=table[!, species_col], kwargs...)
end
