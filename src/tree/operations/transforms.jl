@inline function _tip_edge_index(tree::CompactTree)
    edge_of_child = fill(Int32(0), tree.nnodes)
    @inbounds for edge in 1:tree.nedges
        edge_of_child[Int(tree.child_of_edge[edge])] = Int32(edge)
    end
    return edge_of_child
end

function _compact_tree_with_lengths(tree::CompactTree, edge_length::Vector{Float64}, dist_from_root::Vector{Float64})
    return CompactTree(
        tree.ntips,
        tree.nnodes,
        tree.nedges,
        copy(tree.parent_of_edge),
        copy(tree.child_of_edge),
        edge_length,
        tree.root,
        copy(tree.parent_of_node),
        dist_from_root,
        copy(tree.is_tip),
        copy(tree.tip_ids),
        copy(tree.postorder),
        copy(tree.postorder_internal),
        copy(tree.preorder),
        [copy(ch) for ch in tree.children],
        copy(tree.first_child_edge),
        copy(tree.last_child_edge),
        copy(tree.tipname_to_id),
        copy(tree.tip_labels),
        copy(tree.node_labels),
    )
end

function make_ultrametric(tree::CompactTree)
    tree.ntips > 0 || throw(ArgumentError("tree must contain at least one tip"))
    target_depth = maximum(tree.dist_from_root[tree.tip_ids])
    tip_depths = tree.dist_from_root[tree.tip_ids]
    all(isfinite, tip_depths) || throw(ArgumentError("tree contains non-finite tip depths"))

    edge_length = copy(tree.edge_length)
    dist_from_root = copy(tree.dist_from_root)
    edge_of_child = _tip_edge_index(tree)

    @inbounds for tip in tree.tip_ids
        child = Int(tip)
        edge = Int(edge_of_child[child])
        edge == 0 && throw(ArgumentError("failed to recover terminal edge for tip node $child"))
        add_length = target_depth - dist_from_root[child]
        add_length < 0.0 && throw(ArgumentError("tree tip depth exceeds target depth"))
        edge_length[edge] += add_length
        dist_from_root[child] = target_depth
    end

    return _compact_tree_with_lengths(tree, edge_length, dist_from_root)
end

function make_ultrametric(tree::T) where {T <: NewickTree.Node}
    return from_compact_tree(make_ultrametric(compact_tree(tree)))
end
