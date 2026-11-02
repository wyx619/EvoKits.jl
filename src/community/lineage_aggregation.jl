"""
    LineageCommunityResult

`community` is a CSC `CommunityData` whose columns are the groups recorded in
the owned `grouping::CladeGroupingResult`. `statistic` is :presence, :richness
or :abundance; `normalize` records whether row proportions were requested.
"""
struct LineageCommunityResult{D<:CommunityData}
    community::D
    grouping::CladeGroupingResult
    statistic::Symbol
    normalize::Bool
end

function _copy_clade_grouping(grouping)
    return CladeGroupingResult(copy(grouping.tip_labels), copy(grouping.membership),
        copy(grouping.group_labels), copy(grouping.group_nodes), copy(grouping.group_sizes),
        grouping.cut, grouping.boundary, grouping.tree_height, grouping.requested_k, grouping.actual_k)
end

function _lineage_mapping(grouping)
    labels = _community_data_labels(grouping.tip_labels, length(grouping.membership), "species")
    groups = _community_data_labels(grouping.group_labels, grouping.actual_k, "group")
    length(grouping.group_nodes) == length(grouping.group_sizes) == grouping.actual_k ||
        throw(ArgumentError("grouping metadata dimensions do not match"))
    counts = zeros(Int, grouping.actual_k)
    mapping = Dict{String, String}()
    for (label, group) in zip(labels, grouping.membership)
        1 <= group <= grouping.actual_k || throw(ArgumentError("grouping membership is out of range"))
        mapping[label] = groups[group]
        counts[group] += 1
    end
    all(>(0), counts) && counts == grouping.group_sizes || throw(ArgumentError("grouping sizes do not match membership"))
    return mapping
end

function _collapse_lineages(matrix, grouping, statistic, normalize, species_labels, sample_labels)
    statistic in (:presence, :richness, :abundance) ||
        throw(ArgumentError("statistic must be :presence, :richness or :abundance"))
    species, samples = _community_aggregation_input(matrix, species_labels, sample_labels)
    mapping = _lineage_mapping(grouping)
    unknown = filter(label -> !haskey(mapping, label), species)
    isempty(unknown) || throw(ArgumentError("community species absent from grouping tree: $(join(first(unknown, 5), ", "))"))
    data = aggregate_species(matrix, mapping; group_labels = grouping.group_labels,
        species_labels = species, sample_labels = samples,
        abundance_weighted = statistic === :abundance)
    statistic === :presence && fill!(nonzeros(data.matrix), 1)
    if normalize
        matrix = _community_normalize_groups(data.matrix)
        data = CommunityData(matrix, data.sample_labels, data.species_labels, Val(:validated))
    end
    return LineageCommunityResult(data, _copy_clade_grouping(grouping), statistic, normalize)
end

"""
    collapse_range(data::CommunityData, tree::CompactTree; cut=nothing, k=nothing,
        boundary=nothing, atol=0.0, rtol=1e-8, statistic=:presence, normalize=false)
    collapse_range(data::CommunityData, grouping::CladeGroupingResult;
        statistic=:presence, normalize=false)
    collapse_range(matrix, tree; tip_labels=tree.tip_labels, sample_labels=nothing, ...)
    collapse_range(matrix, grouping; tip_labels=grouping.tip_labels, sample_labels=nothing, ...)

Collapse species columns into tree-derived lineage groups. A tree input calls
`get_clades` once; pass its result to reuse the same full-tree mapping across
community inputs. Return `LineageCommunityResult` with labeled CSC `community`
and an independently owned `grouping` containing the complete tip-to-group map.

Default `statistic=:presence` counts each lineage at most once per site (0/1).
`:richness` counts its present species, and `:abundance` sums original values,
preserving their type (Bool becomes Int), with overflow checks. `normalize=true`
returns per-site Float64 proportions *after* the selected statistic. Empty
sites stay zero. Pass `community` to `beta_diss`, then its result to `phyloregion`.

Species are matched exactly by name. Missing tree tips and their groups remain
as zero contributions/columns; every group of the complete tree is retained.
Any unknown species is an error, even in an all-zero column. No tree pruning,
silent intersection, name rewriting, dense padding, dense membership matrix,
or long occurrence table is used. Matrix labels default to the tree/grouping
tip order and therefore require one column per tree tip unless specified.
Input matrices, labels and grouping arrays are not mutated.
"""
function collapse_range(data::CommunityData, grouping::CladeGroupingResult;
    statistic::Symbol = :presence, normalize::Bool = false,
)
    _validate_community_data(data)
    return _collapse_lineages(data.matrix, grouping, statistic, normalize,
        data.species_labels, data.sample_labels)
end

function collapse_range(matrix::AbstractMatrix{<:Real}, grouping::CladeGroupingResult;
    statistic::Symbol = :presence, normalize::Bool = false,
    tip_labels = grouping.tip_labels, sample_labels = nothing,
)
    return _collapse_lineages(matrix, grouping, statistic, normalize, tip_labels, sample_labels)
end

function collapse_range(data::CommunityData, tree::CompactTree; cut = nothing,
    k = nothing, boundary = nothing, atol = 0.0, rtol = 1e-8,
    statistic::Symbol = :presence, normalize::Bool = false,
)
    grouping = get_clades(tree; cut, k, boundary, atol, rtol)
    return collapse_range(data, grouping; statistic, normalize)
end

function collapse_range(matrix::AbstractMatrix{<:Real}, tree::CompactTree; cut = nothing,
    k = nothing, boundary = nothing, atol = 0.0, rtol = 1e-8,
    statistic::Symbol = :presence, normalize::Bool = false,
    tip_labels = tree.tip_labels, sample_labels = nothing,
)
    grouping = get_clades(tree; cut, k, boundary, atol, rtol)
    return collapse_range(matrix, grouping; statistic, normalize, tip_labels, sample_labels)
end
