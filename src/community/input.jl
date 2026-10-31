@inline function _tools_sample_labels(n::Integer, labels::Union{Nothing, AbstractVector{<:AbstractString}})
    labels === nothing && return ["sample_$i" for i in 1:n]
    length(labels) == n || throw(ArgumentError("sample_labels length must match number of sample rows"))
    return String.(labels)
end

function _tools_tip_permutation(tree::CompactTree, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}})
    tip_labels === nothing && return collect(1:tree.ntips)
    length(tip_labels) == tree.ntips || throw(ArgumentError("tip_labels length must match number of sample columns"))
    labels = String.(tip_labels)
    length(unique(labels)) == length(labels) || throw(ArgumentError("tip_labels must not contain duplicates"))
    tip_position = zeros(Int32, tree.nnodes)
    for (position, node) in enumerate(tree.tip_ids)
        tip_position[node] = Int32(position)
    end
    perm = Vector{Int}(undef, tree.ntips)
    for (j, label) in enumerate(labels)
        node = get(tree.tipname_to_id, label, Int32(0))
        node == 0 && throw(ArgumentError("sample column label $label is not present in tree tips"))
        idx = tip_position[node]
        idx == 0 && throw(ArgumentError("tree tip index not found for label $label"))
        perm[j] = idx
    end
    return perm
end

function _tools_tip_permutation(tree::CompactTree, samp::AbstractDataFrame, tip_labels::Nothing)
    labels = String.(names(samp))
    length(labels) == tree.ntips || throw(ArgumentError("dataframe must have one column per tree tip"))
    return _tools_tip_permutation(tree, labels)
end

function _tools_tip_permutation(tree::CompactTree, samp::AbstractDataFrame, tip_labels::AbstractVector{<:AbstractString})
    length(tip_labels) == ncol(samp) || throw(ArgumentError("tip_labels length must match number of sample columns"))
    return _tools_tip_permutation(tree, tip_labels)
end

function _tools_abundance_only(samp::AbstractMatrix{<:Real}, perm::Vector{Int})
    all(x -> isfinite(x) && x >= 0, samp) || throw(ArgumentError("community values must be finite and non-negative"))
    abund = Matrix{Float64}(undef, size(samp, 1), length(perm))
    @inbounds for i in axes(samp, 1), j in axes(samp, 2)
        abund[i, perm[j]] = Float64(samp[i, j])
    end
    return abund
end

function _tools_presence_only(samp::AbstractMatrix{<:Real}, perm::Vector{Int})
    all(x -> isfinite(x) && x >= 0, samp) || throw(ArgumentError("community values must be finite and non-negative"))
    presence = Matrix{Bool}(undef, size(samp, 1), length(perm))
    @inbounds for i in axes(samp, 1), j in axes(samp, 2)
        presence[i, perm[j]] = samp[i, j] > 0
    end
    return presence
end

function _tools_presence_matrix(samp::AbstractMatrix{<:Real}, perm::Vector{Int})
    return _tools_presence_only(samp, perm)
end

function _community_abundance(samp::AbstractMatrix{<:Real}, tree::CompactTree, tip_labels)
    size(samp, 2) == tree.ntips || throw(ArgumentError("community matrix must have one column per tree tip"))
    perm = _tools_tip_permutation(tree, tip_labels)
    return _tools_abundance_only(samp, perm)
end

function _community_presence(samp::AbstractMatrix{<:Real}, tree::CompactTree, tip_labels)
    abundance = _community_abundance(samp, tree, tip_labels)
    return abundance .> 0
end

# The presence-based kernels below visit stored entries of CSC input directly.
# Legacy SES helpers above intentionally keep their dense randomization buffers.
function _validate_community_values(community::AbstractMatrix{<:Real})
    Base.require_one_based_indexing(community)
    values = community isa SparseMatrixCSC ? nonzeros(community) : community
    all(x -> isfinite(x) && x >= 0, values) ||
        throw(ArgumentError("community values must be finite and non-negative"))
    return nothing
end

function _validate_community_branch_lengths(tree::CompactTree)
    all(x -> isfinite(x) && x >= 0, tree.edge_length) ||
        throw(ArgumentError("community metrics require finite, non-negative branch lengths"))
    isfinite(sum(tree.edge_length)) ||
        throw(ArgumentError("total branch length must be finite"))
    return nothing
end

function _community_presence_rows(community::AbstractMatrix{<:Real}, perm::AbstractVector{<:Integer})
    _validate_community_values(community)
    rows = [Int32[] for _ in axes(community, 1)]
    for col in axes(community, 2), site in axes(community, 1)
        community[site, col] > 0 && push!(rows[site], Int32(perm[col]))
    end
    return rows
end

function _community_presence_rows(community::SparseMatrixCSC{<:Real}, perm::AbstractVector{<:Integer})
    _validate_community_values(community)
    rows = [Int32[] for _ in axes(community, 1)]
    row_indices, values = rowvals(community), nonzeros(community)
    for col in axes(community, 2), entry in nzrange(community, col)
        values[entry] > 0 && push!(rows[row_indices[entry]], Int32(perm[col]))
    end
    return rows
end

function _community_tip_rows(community::AbstractMatrix{<:Real}, tree::CompactTree, tip_labels)
    size(community, 2) == tree.ntips ||
        throw(ArgumentError("community matrix must have one column per tree tip"))
    perm = _tools_tip_permutation(tree, tip_labels)
    rows = _community_presence_rows(community, perm)
    # Match tree-tip order even when input columns were permuted.
    if !issorted(perm)
        foreach(sort!, rows)
    end
    return rows
end

function _community_edge_index(tree::CompactTree)
    edge_of_child = zeros(Int32, tree.nnodes)
    for edge in 1:tree.nedges
        edge_of_child[tree.child_of_edge[edge]] = Int32(edge)
    end
    return edge_of_child
end

# A branch is visited once per site. Stop at the first visited ancestor, so
# neither a dense site-by-edge matrix nor all tip-to-root paths are needed.
function _walk_community_branches(f, rows, tree::CompactTree;
    roots::Union{Nothing, Vector{Int32}} = nothing,
    edge_of_child::Vector{Int32} = _community_edge_index(tree),
)
    seen = zeros(Int, tree.nnodes)
    @inbounds for (site, tips) in enumerate(rows)
        root = roots === nothing ? tree.root : roots[site]
        for tip in tips
            node = tree.tip_ids[tip]
            while node != root && seen[node] != site
                seen[node] = site
                edge = edge_of_child[node]
                f(site, edge)
                node = tree.parent_of_node[node]
            end
        end
    end
    return nothing
end

function _community_branch_sets(rows, tree::CompactTree)
    branches = [Int32[] for _ in eachindex(rows)]
    _walk_community_branches(rows, tree) do site, edge
        push!(branches[site], edge)
    end
    foreach(sort!, branches)
    return branches
end
