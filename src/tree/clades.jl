function descendant_tips(tree::CompactTree, node; labels::Bool = false, by::Symbol = :node)
    node_id = _tree_resolve_node(tree, node; by = by)
    positions = Int[]
    @inbounds for (position, tip) in enumerate(tree.tip_ids)
        current = Int(tip)
        while current != Int(tree.root) && current != Int(node_id)
            current = Int(tree.parent_of_node[current])
        end
        current == Int(node_id) && push!(positions, position)
    end
    return labels ? tree.tip_labels[positions] : positions
end

function clade_members(tree::CompactTree; include_nodes::Bool = false)
    members = [Int[] for _ in 1:tree.nnodes]
    @inbounds for node32 in tree.postorder
        node = Int(node32)
        if tree.is_tip[node]
            position = findfirst(==(node32), tree.tip_ids)
            position === nothing || push!(members[node], position)
        else
            for child in tree.children[node]
                append!(members[node], members[Int(child)])
            end
        end
    end
    if include_nodes
        return members
    end
    return [members[Int(node)] for node in tree.postorder_internal]
end

function subtree_at_node(tree::CompactTree, node; by::Symbol = :node)
    labels = descendant_tips(tree, node; labels = true, by = by)
    return keep_tip(tree, labels)
end
