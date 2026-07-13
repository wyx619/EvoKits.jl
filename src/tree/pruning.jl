function _normalize_tip_set(tree::CompactTree, tips)
    labels = tips isa AbstractString ? [String(tips)] : String.(collect(tips))
    length(unique(labels)) == length(labels) || throw(ArgumentError("tips must not contain duplicates"))
    for label in labels
        haskey(tree.tipname_to_id, label) || throw(ArgumentError("tip $label is not present in the tree"))
    end
    return Set(labels)
end

function _edge_lengths_to_node(tree::CompactTree)
    lengths = zeros(Float64, tree.nnodes)
    @inbounds for edge in 1:tree.nedges
        lengths[Int(tree.child_of_edge[edge])] = tree.edge_length[edge]
    end
    return lengths
end

function keep_tip(tree::CompactTree, tips)
    keep = _normalize_tip_set(tree, tips)
    isempty(keep) && throw(ArgumentError("Cannot keep zero tips from a tree"))
    length(keep) == tree.ntips && return tree

    retained_descendants = falses(tree.nnodes)
    @inbounds for node in tree.postorder
        inode = Int(node)
        if tree.is_tip[inode]
            retained_descendants[inode] = tree.node_labels[inode] in keep
        else
            retained_descendants[inode] = any(retained_descendants[Int(child)] for child in tree.children[inode])
        end
    end

    edge_lengths = _edge_lengths_to_node(tree)
    function collapse_unary(node::Int, length_from_parent::Float64)
        retained_children = Int[child for child in tree.children[node] if retained_descendants[Int(child)]]
        if tree.is_tip[node] || length(retained_children) != 1
            return node, length_from_parent
        end
        child = only(retained_children)
        return collapse_unary(child, length_from_parent + edge_lengths[child])
    end

    root, _ = collapse_unary(Int(tree.root), 0.0)
    parent_of_node = Int32[]
    children = Vector{Vector{Int32}}()
    edge_length_to_node = Float64[]
    is_tip = BitVector()
    tip_labels = String[]
    node_labels = String[]

    function add_node(old_node::Int, parent::Int32, length_from_parent::Float64)
        new_node = Int32(length(parent_of_node) + 1)
        push!(parent_of_node, parent)
        push!(children, Int32[])
        push!(edge_length_to_node, length_from_parent)
        push!(is_tip, tree.is_tip[old_node])
        push!(node_labels, tree.node_labels[old_node])
        parent != 0 && push!(children[Int(parent)], new_node)
        if tree.is_tip[old_node]
            push!(tip_labels, tree.node_labels[old_node])
            return new_node
        end
        for child0 in tree.children[old_node]
            child = Int(child0)
            retained_descendants[child] || continue
            next_node, next_length = collapse_unary(child, edge_lengths[child])
            add_node(next_node, new_node, next_length)
        end
        return new_node
    end

    add_node(root, Int32(0), 0.0)
    return _compact_tree_from_topology(
        parent_of_node,
        children,
        edge_length_to_node,
        is_tip,
        tip_labels,
        node_labels,
    )
end

function drop_tip(tree::CompactTree, tips)
    drop = _normalize_tip_set(tree, tips)
    keep = setdiff(Set(tree.tip_labels), drop)
    isempty(keep) && throw(ArgumentError("Cannot drop all tips from a tree"))
    return keep_tip(tree, keep)
end
