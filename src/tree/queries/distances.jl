function _tree_resolve_node(tree::CompactTree, value; by::Symbol = :tip)
    if value isa AbstractString
        node = get(tree.tipname_to_id, String(value), Int32(0))
        node == 0 && throw(ArgumentError("unknown tree tip: $value"))
        return node
    elseif by === :tip
        1 <= value <= tree.ntips || throw(ArgumentError("tip position $value is outside 1:$(tree.ntips)"))
        return tree.tip_ids[Int(value)]
    elseif by === :node
        1 <= value <= tree.nnodes || throw(ArgumentError("node id $value is outside 1:$(tree.nnodes)"))
        return Int32(value)
    else
        throw(ArgumentError("by must be :tip or :node"))
    end
end

function mrca(tree::CompactTree, first, second; by::Symbol = :tip)
    index = _build_tree_query_index(tree)
    return Int(_tree_lca(index, _tree_resolve_node(tree, first; by = by), _tree_resolve_node(tree, second; by = by)))
end

function mrca(tree::CompactTree, values::AbstractVector; by::Symbol = :tip)
    isempty(values) && throw(ArgumentError("mrca requires at least one tip or node"))
    nodes = [_tree_resolve_node(tree, value; by = by) for value in values]
    index = _build_tree_query_index(tree)
    result = first(nodes)
    for node in Iterators.drop(nodes, 1)
        result = _tree_lca(index, result, node)
    end
    return Int(result)
end

function distance(tree::CompactTree, first, second; weighted::Bool = true, by::Symbol = :tip)
    first_node = _tree_resolve_node(tree, first; by = by)
    second_node = _tree_resolve_node(tree, second; by = by)
    index = _build_tree_query_index(tree)
    ancestor = _tree_lca(index, first_node, second_node)
    if weighted
        return tree.dist_from_root[Int(first_node)] + tree.dist_from_root[Int(second_node)] - 2.0 * tree.dist_from_root[Int(ancestor)]
    end
    depth = zeros(Int, tree.nnodes)
    @inbounds for node in tree.preorder
        inode = Int(node)
        parent = Int(tree.parent_of_node[inode])
        depth[inode] = parent == 0 ? 0 : depth[parent] + 1
    end
    return depth[Int(first_node)] + depth[Int(second_node)] - 2 * depth[Int(ancestor)]
end

function pairwise_distance(tree::CompactTree; tips = nothing, weighted::Bool = true, by::Symbol = :tip)
    values = tips === nothing ? collect(1:tree.ntips) : collect(tips)
    isempty(values) && return zeros(Float64, 0, 0)
    n = length(values)
    nodes = [_tree_resolve_node(tree, value; by = by) for value in values]
    index = _build_tree_query_index(tree)
    result = zeros(Float64, n, n)
    depth = if weighted
        nothing
    else
        values = zeros(Int, tree.nnodes)
        @inbounds for node in tree.preorder
            inode = Int(node)
            parent = Int(tree.parent_of_node[inode])
            values[inode] = parent == 0 ? 0 : values[parent] + 1
        end
        values
    end
    @inbounds for i in 1:n
        for j in 1:(i - 1)
            ancestor = _tree_lca(index, nodes[i], nodes[j])
            result[i, j] = result[j, i] = if weighted
                tree.dist_from_root[Int(nodes[i])] + tree.dist_from_root[Int(nodes[j])] - 2.0 * tree.dist_from_root[Int(ancestor)]
            else
                depth[Int(nodes[i])] + depth[Int(nodes[j])] - 2 * depth[Int(ancestor)]
            end
        end
    end
    return result
end

function distances_from_root(tree::CompactTree; weighted::Bool = true)
    weighted && return copy(tree.dist_from_root)
    depth = zeros(Int, tree.nnodes)
    @inbounds for node in tree.preorder
        inode = Int(node)
        parent = Int(tree.parent_of_node[inode])
        depth[inode] = parent == 0 ? 0 : depth[parent] + 1
    end
    return depth
end

function distances_to_tip(tree::CompactTree, focal; weighted::Bool = true, by::Symbol = :tip)
    focal_node = _tree_resolve_node(tree, focal; by = by)
    nodes = tree.tip_ids
    index = _build_tree_query_index(tree)
    result = Vector{Float64}(undef, tree.ntips)
    depth = weighted ? nothing : distances_from_root(tree; weighted = false)
    @inbounds for (i, node) in enumerate(nodes)
        ancestor = _tree_lca(index, focal_node, node)
        result[i] = weighted ?
            tree.dist_from_root[Int(focal_node)] + tree.dist_from_root[Int(node)] - 2.0 * tree.dist_from_root[Int(ancestor)] :
            Float64(depth[Int(focal_node)] + depth[Int(node)] - 2 * depth[Int(ancestor)])
    end
    return result
end

function _tree_shared_path_covariance(tree::CompactTree)
    index = _build_tree_query_index(tree)
    covariance = Matrix{Float64}(undef, tree.ntips, tree.ntips)
    @inbounds for i in 1:tree.ntips
        for j in 1:i
            ancestor = _tree_lca(index, tree.tip_ids[i], tree.tip_ids[j])
            value = tree.dist_from_root[Int(ancestor)]
            covariance[i, j] = value
            covariance[j, i] = value
        end
    end
    return covariance
end
