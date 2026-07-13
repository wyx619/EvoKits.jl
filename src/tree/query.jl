struct _TreeQueryIndex
    tin::Vector{Int32}
    tout::Vector{Int32}
    up::Matrix{Int32}
end

function _build_tree_query_index(tree::CompactTree)
    nnodes = tree.nnodes
    tin = zeros(Int32, nnodes)
    subtree_size = ones(Int32, nnodes)
    @inbounds for (rank, node) in enumerate(tree.preorder)
        tin[Int(node)] = Int32(rank)
    end
    @inbounds for node in tree.postorder
        inode = Int(node)
        tree.is_tip[inode] && continue
        subtree_size[inode] = Int32(1 + sum(subtree_size[Int(child)] for child in tree.children[inode]))
    end
    tout = tin .+ subtree_size .- Int32(1)

    levels = max(1, ceil(Int, log2(max(nnodes, 2))))
    up = Matrix{Int32}(undef, levels, nnodes)
    root = Int(tree.root)
    @inbounds for node in 1:nnodes
        parent = Int(tree.parent_of_node[node])
        up[1, node] = Int32(parent == 0 ? root : parent)
    end
    @inbounds for level in 2:levels, node in 1:nnodes
        up[level, node] = up[level - 1, Int(up[level - 1, node])]
    end
    return _TreeQueryIndex(tin, tout, up)
end

@inline function _tree_is_ancestor(index::_TreeQueryIndex, ancestor::Int32, node::Int32)
    return index.tin[Int(ancestor)] <= index.tin[Int(node)] && index.tout[Int(node)] <= index.tout[Int(ancestor)]
end

function _tree_lca(index::_TreeQueryIndex, first::Int32, second::Int32)
    _tree_is_ancestor(index, first, second) && return first
    _tree_is_ancestor(index, second, first) && return second
    node = first
    @inbounds for level in size(index.up, 1):-1:1
        candidate = index.up[level, Int(node)]
        !_tree_is_ancestor(index, candidate, second) && (node = candidate)
    end
    return index.up[1, Int(node)]
end

function _community_tip_positions(tree::CompactTree, tip_positions::AbstractVector{<:Integer})
    positions = Int.(tip_positions)
    length(positions) >= 2 || throw(ArgumentError("a community must contain at least two tips"))
    all((1 .<= positions) .& (positions .<= tree.ntips)) || throw(ArgumentError("community tip positions are outside 1:$(tree.ntips)"))
    length(unique(positions)) == length(positions) || throw(ArgumentError("community contains duplicate tips"))
    return positions
end

function _community_metrics(tree::CompactTree, tip_positions::AbstractVector{<:Integer}, index::_TreeQueryIndex)
    positions = _community_tip_positions(tree, tip_positions)
    selected = Int32[tree.tip_ids[position] for position in positions]
    sort!(selected; by = node -> index.tin[Int(node)])
    virtual_nodes = copy(selected)
    for i in 1:(length(selected) - 1)
        push!(virtual_nodes, _tree_lca(index, selected[i], selected[i + 1]))
    end
    sort!(virtual_nodes; by = node -> index.tin[Int(node)])
    unique!(virtual_nodes)

    virtual_index = Dict{Int32, Int32}(node => Int32(i) for (i, node) in enumerate(virtual_nodes))
    marked = falses(length(virtual_nodes))
    @inbounds for node in selected
        marked[Int(virtual_index[node])] = true
    end

    parent_index = zeros(Int32, length(virtual_nodes))
    stack = Int32[]
    @inbounds for i in eachindex(virtual_nodes)
        node = virtual_nodes[i]
        while !isempty(stack) && !_tree_is_ancestor(index, virtual_nodes[Int(stack[end])], node)
            pop!(stack)
        end
        !isempty(stack) && (parent_index[i] = stack[end])
        push!(stack, Int32(i))
    end
    children = [Int32[] for _ in eachindex(virtual_nodes)]
    @inbounds for i in 2:length(virtual_nodes)
        push!(children[Int(parent_index[i])], Int32(i))
    end

    richness = length(selected)
    counts = Int32.(marked)
    down = fill(Inf, length(virtual_nodes))
    down[marked] .= 0.0
    pair_sum = 0.0
    @inbounds for i in length(virtual_nodes):-1:2
        parent = Int(parent_index[i])
        edge_length = tree.dist_from_root[Int(virtual_nodes[i])] - tree.dist_from_root[Int(virtual_nodes[parent])]
        pair_sum += edge_length * counts[i] * (richness - counts[i])
        counts[parent] += counts[i]
        down[parent] = min(down[parent], down[i] + edge_length)
    end

    updist = fill(Inf, length(virtual_nodes))
    @inbounds for parent in eachindex(virtual_nodes)
        best1 = marked[parent] ? 0.0 : Inf
        best2 = Inf
        best_child = Int32(0)
        for child32 in children[parent]
            child = Int(child32)
            edge_length = tree.dist_from_root[Int(virtual_nodes[child])] - tree.dist_from_root[Int(virtual_nodes[parent])]
            candidate = down[child] + edge_length
            if candidate < best1
                best2, best1, best_child = best1, candidate, child32
            elseif candidate < best2
                best2 = candidate
            end
        end
        for child32 in children[parent]
            child = Int(child32)
            edge_length = tree.dist_from_root[Int(virtual_nodes[child])] - tree.dist_from_root[Int(virtual_nodes[parent])]
            sibling_best = child32 == best_child ? best2 : best1
            updist[child] = edge_length + min(updist[parent], sibling_best)
        end
    end

    nearest_sum = 0.0
    @inbounds for i in eachindex(virtual_nodes)
        marked[i] || continue
        child_best = Inf
        for child32 in children[i]
            child = Int(child32)
            edge_length = tree.dist_from_root[Int(virtual_nodes[child])] - tree.dist_from_root[Int(virtual_nodes[i])]
            child_best = min(child_best, down[child] + edge_length)
        end
        nearest_sum += min(updist[i], child_best)
    end

    mpd = pair_sum / (richness * (richness - 1) / 2)
    mntd = nearest_sum / richness
    return mpd, mntd
end

function _community_metrics(tree::CompactTree, tip_positions::AbstractVector{<:Integer})
    return _community_metrics(tree, tip_positions, _build_tree_query_index(tree))
end

function _tools_tip_paths(tree::CompactTree)
    edge_of_child = zeros(Int32, tree.nnodes)
    @inbounds for edge in 1:tree.nedges
        edge_of_child[Int(tree.child_of_edge[edge])] = Int32(edge)
    end
    paths = Vector{Vector{Int32}}(undef, tree.ntips)
    @inbounds for (i, tip) in enumerate(tree.tip_ids)
        node = Int(tip)
        path = Int32[]
        while node != Int(tree.root)
            edge = edge_of_child[node]
            edge == 0 && break
            push!(path, edge)
            node = Int(tree.parent_of_node[node])
        end
        paths[i] = path
    end
    return paths
end

function _tools_descendant_tip_counts(tree::CompactTree)
    counts = zeros(Int32, tree.nnodes)
    @inbounds for node in tree.postorder
        inode = Int(node)
        if tree.is_tip[inode]
            counts[inode] = 1
        else
            total = Int32(0)
            for child in tree.children[inode]
                total += counts[Int(child)]
            end
            counts[inode] = total
        end
    end
    return counts
end

@inline function _tools_mrca_depth_from_paths(paths::Vector{Vector{Int32}}, present::Vector{Int}, tree::CompactTree)
    length(present) <= 1 && return 0.0
    ref = paths[present[1]]
    shared = 0
    max_shared = length(ref)
    for pos in 1:max_shared
        edge = ref[end - pos + 1]
        ok = true
        for k in 2:length(present)
            path = paths[present[k]]
            if pos > length(path) || path[end - pos + 1] != edge
                ok = false
                break
            end
        end
        ok || break
        shared = pos
    end
    shared == 0 && return 0.0
    mrca_edge = ref[end - shared + 1]
    return tree.dist_from_root[Int(tree.child_of_edge[Int(mrca_edge)])]
end
