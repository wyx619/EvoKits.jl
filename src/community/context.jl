"""
    CommunityContext

Reusable presence-based community snapshot created by [`prepare_community`](@ref).
Use it with `pd`, `pd_dataframe`, `weighted_endemism`, `phylo_endemism`,
`relative_pd`, `phylobeta`, `unifrac`, and `phylosor`. The context owns its
preprocessed arrays; its fields are internal and must not be mutated.
"""
struct CommunityContext
    sample_labels::Vector{String}
    branches::Vector{Vector{Int32}}
    branch_lengths::Vector{Float64}
    richness::Vector{Int}
    species_occupancy::Vector{Int}
    branch_occupancy::Vector{Int}
    pd_values::Vector{Float64}
    pd_no_root::Vector{Float64}
    positive_counts::Vector{Int}
    we_values::Vector{Float64}
    pe_weighted::Vector{Float64}
    pe_strict::Vector{Float64}
    total_branch_length::Float64
    positive_branch_count::Int
end

"""
    prepare_community(community, tree; tip_labels=nothing, sample_labels=nothing)

Validate and snapshot a complete community collection and its tree for repeated
presence-based calculations. Accepts dense matrices and `SparseMatrixCSC`;
rows are sites, columns are species, and finite positive values denote presence.
Supply one column per tree tip, including absent species as zero columns.

Returns a `CommunityContext` that reuses label matching, sorted branch sets,
species/branch occupancy, and per-site summaries. Occupancy always refers to the
complete collection supplied here; `pairs` indices in later `phylobeta` calls
refer to these original rows. Species occupancy retains the input column order.

Later changes to the input matrix, tree, or labels do not affect the context.
Rebuild it to analyze changed data. Metric results own their arrays and may be
modified independently. The context retains neither the input matrix nor the
tree, and never caches full site-by-site matrices. Its storage is proportional
to the per-site branch sets plus site, species, and edge summaries.
"""
function prepare_community(community::AbstractMatrix{<:Real}, tree::CompactTree;
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    _validate_community_branch_lengths(tree)
    size(community, 2) == tree.ntips ||
        throw(ArgumentError("community matrix must have one column per tree tip"))
    labels = _tools_sample_labels(size(community, 1), sample_labels)
    perm = _tools_tip_permutation(tree, tip_labels)
    rows = _community_tip_rows(community, perm)
    lengths = copy(tree.edge_length)
    richness = length.(rows)
    tip_occupancy = zeros(Int, tree.ntips)
    for tips in rows, tip in tips
        tip_occupancy[tip] += 1
    end
    we = zeros(Float64, length(rows))
    for (site, tips) in enumerate(rows), tip in tips
        we[site] += 1 / tip_occupancy[tip]
    end

    branches = [Int32[] for _ in eachindex(rows)]
    occupancy = zeros(Int, tree.nedges)
    pd_values = zeros(Float64, length(rows))
    positive_counts = zeros(Int, length(rows))
    edge_index = _community_edge_index(tree)
    _walk_community_branches(rows, tree; edge_of_child = edge_index) do site, edge
        push!(branches[site], edge)
        occupancy[edge] += 1
        pd_values[site] += lengths[edge]
        positive_counts[site] += lengths[edge] > 0
    end
    foreach(sort!, branches)
    weighted, strict = zeros(Float64, length(rows)), zeros(Float64, length(rows))
    for (site, edges) in enumerate(branches), edge in edges
        weighted[site] += lengths[edge] / occupancy[edge]
        occupancy[edge] == 1 && (strict[site] += lengths[edge])
    end

    # Sum below each MRCA directly instead of subtracting a potentially long
    # shared root path; this preserves short branches near floating cancellation.
    no_root = zeros(Float64, length(rows))
    roots = _community_roots(rows, tree)
    _walk_community_branches(rows, tree; roots = roots, edge_of_child = edge_index) do site, edge
        no_root[site] += lengths[edge]
    end
    no_root[richness .== 1] .= NaN
    return CommunityContext(labels, branches, lengths, richness, tip_occupancy[perm],
        occupancy, pd_values, no_root, positive_counts, we, weighted, strict,
        sum(lengths), count(>(0), lengths))
end

"""`pd(context; include_root=true)` reuses prepared PD and richness summaries."""
function pd(context::CommunityContext; include_root::Bool = true)
    values = include_root ? context.pd_values : context.pd_no_root
    return PDResult(sample_labels = copy(context.sample_labels), pd = copy(values),
        richness = copy(context.richness), include_root = include_root)
end

function pd_dataframe(context::CommunityContext; include_root::Bool = true)
    result = pd(context; include_root = include_root)
    return DataFrame(sample = result.sample_labels, PD = result.pd, SR = result.richness)
end

"""`weighted_endemism(context)` uses occupancy from the complete prepared collection."""
function weighted_endemism(context::CommunityContext)
    return WeightedEndemismResult(copy(context.sample_labels), copy(context.we_values),
        copy(context.richness), copy(context.species_occupancy))
end

"""`phylo_endemism(context; weighted=true)` reuses prepared global branch occupancy."""
function phylo_endemism(context::CommunityContext; weighted::Bool = true)
    values = weighted ? context.pe_weighted : context.pe_strict
    return PhyloEndemismResult(copy(context.sample_labels), copy(values),
        copy(context.richness), copy(context.branch_occupancy), weighted)
end

"""`relative_pd(context)` reuses original PD and equal-positive-branch summaries."""
function relative_pd(context::CommunityContext)
    if context.total_branch_length > 0
        original = context.pd_values ./ context.total_branch_length
        equal = context.positive_counts ./ context.positive_branch_count
    else
        original, equal = (fill(NaN, length(context.richness)) for _ in 1:2)
    end
    rpd = [equal[i] > 0 ? original[i] / equal[i] : NaN for i in eachindex(equal)]
    return RelativePDResult(copy(context.sample_labels), rpd, original, equal,
        copy(context.richness), context.total_branch_length, context.positive_branch_count)
end

"""`phylobeta(context; index_family=:sorensen, pairs=nothing)` reuses prepared branch sets."""
function phylobeta(context::CommunityContext; index_family::Symbol = :sorensen, pairs = nothing)
    _validate_community_beta_family(index_family)
    selected = _community_pairs(pairs, length(context.richness))
    total, turnover, nestedness = _phylobeta(context.branches, context.branch_lengths,
        index_family, selected)
    return PhylobetaResult(copy(context.sample_labels), total, turnover, nestedness,
        index_family, selected)
end

"""`unifrac(context)` calculates distances using the prepared branch sets."""
unifrac(context::CommunityContext) = _unifrac(context.branches, context.branch_lengths)

"""`phylosor(context)` calculates similarity using the prepared branch sets."""
phylosor(context::CommunityContext) = _phylosor(context.branches, context.branch_lengths)
