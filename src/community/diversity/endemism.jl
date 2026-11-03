"""Weighted endemism by site; occupancy follows the input species-column order."""
struct WeightedEndemismResult
    sample_labels::Vector{String}
    we::Vector{Float64}
    richness::Vector{Int}
    species_occupancy::Vector{Int}
end

"""Phylogenetic endemism by site; occupancy follows the tree's edge order."""
struct PhyloEndemismResult
    sample_labels::Vector{String}
    pe::Vector{Float64}
    richness::Vector{Int}
    branch_occupancy::Vector{Int}
    weighted::Bool
end

"""
    weighted_endemism(community; sample_labels=nothing)

Sum reciprocal species occupancy over each site's present species. Rows are
sites, columns are species, and positive abundance denotes presence. Supports
dense matrices and SparseMatrixCSC without densifying sparse input. Empty
sites return zero. Occupancy is calculated over the complete input collection.
"""
function weighted_endemism(community::AbstractMatrix{<:Real};
    sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    labels = _tools_sample_labels(size(community, 1), sample_labels)
    rows = _community_presence_rows(community, axes(community, 2))
    occupancy = zeros(Int, size(community, 2))
    for row in rows, tip in row
        occupancy[tip] += 1
    end
    values = zeros(Float64, length(rows))
    for (site, row) in enumerate(rows), tip in row
        values[site] += 1 / occupancy[tip]
    end
    return WeightedEndemismResult(labels, values, length.(rows), occupancy)
end

"""
    phylo_endemism(community, tree; weighted=true, tip_labels=nothing, sample_labels=nothing)

Sum branch length divided by branch occupancy (weighted=true), or lengths of
branches restricted to one site (weighted=false). Each site's subtree includes
its paths to the root. Empty sites return zero. All sites must be supplied
together because branch occupancy is measured over the complete collection.
"""
function phylo_endemism(community::AbstractMatrix{<:Real}, tree::CompactTree;
    weighted::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    _validate_community_branch_lengths(tree)
    labels = _tools_sample_labels(size(community, 1), sample_labels)
    tips = _community_tip_rows(community, tree, tip_labels)
    occupancy = zeros(Int, tree.nedges)
    _walk_community_branches(tips, tree) do site, edge
        occupancy[edge] += 1
    end
    values = zeros(Float64, length(tips))
    _walk_community_branches(tips, tree) do site, edge
        if weighted
            values[site] += tree.edge_length[edge] / occupancy[edge]
        elseif occupancy[edge] == 1
            values[site] += tree.edge_length[edge]
        end
    end
    return PhyloEndemismResult(labels, values, length.(tips), occupancy, weighted)
end
