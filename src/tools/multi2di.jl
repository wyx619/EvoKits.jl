Base.@kwdef mutable struct _Multi2DiItem
    root::Int32
    length::Float64
end

function _compact_tree_from_topology(parent_of_node::Vector{Int32}, children::Vector{Vector{Int32}}, edge_length_to_node::Vector{Float64}, is_tip::BitVector, tip_labels::Vector{String}, node_labels::Vector{String})
    nnodes = length(parent_of_node)
    roots = findall(==(Int32(0)), parent_of_node)
    length(roots) == 1 || throw(ArgumentError("tree must contain exactly one root"))
    root = Int32(only(roots))

    dist_from_root = zeros(Float64, nnodes)
    preorder = Int32[]
    stack = Int32[root]
    while !isempty(stack)
        node = pop!(stack)
        push!(preorder, node)
        for child in Iterators.reverse(children[Int(node)])
            dist_from_root[Int(child)] = dist_from_root[Int(node)] + edge_length_to_node[Int(child)]
            push!(stack, child)
        end
    end

    postorder = Int32[]
    postorder_internal = Int32[]
    function dfs(node::Int32)
        for child in children[Int(node)]
            dfs(child)
        end
        push!(postorder, node)
        !is_tip[Int(node)] && push!(postorder_internal, node)
        return nothing
    end
    dfs(root)

    tip_ids = Int32[findall(is_tip)...]
    ntips = length(tip_ids)
    nedges = nnodes - 1
    tip_labels_by_node = fill("", nnodes)
    for (i, tip) in enumerate(tip_ids)
        tip_labels_by_node[Int(tip)] = tip_labels[i]
    end
    ordered_tip_labels = [tip_labels_by_node[Int(tip)] for tip in tip_ids]
    tipname_to_id = Dict{String, Int32}(ordered_tip_labels[i] => tip_ids[i] for i in eachindex(tip_ids))

    parent_of_edge = Vector{Int32}(undef, nedges)
    child_of_edge = Vector{Int32}(undef, nedges)
    edge_length = Vector{Float64}(undef, nedges)
    first_child_edge = fill(Int32(0), nnodes)
    last_child_edge = fill(Int32(0), nnodes)

    edge_idx = Int32(1)
    for parent in Int32.(1:nnodes)
        for child in children[Int(parent)]
            parent_of_edge[edge_idx] = parent
            child_of_edge[edge_idx] = child
            edge_length[edge_idx] = edge_length_to_node[Int(child)]
            if first_child_edge[Int(parent)] == 0
                first_child_edge[Int(parent)] = edge_idx
            end
            last_child_edge[Int(parent)] = edge_idx
            edge_idx += 1
        end
    end

    return CompactTree(
        ntips,
        nnodes,
        nedges,
        parent_of_edge,
        child_of_edge,
        edge_length,
        root,
        parent_of_node,
        dist_from_root,
        is_tip,
        tip_ids,
        postorder,
        postorder_internal,
        preorder,
        children,
        first_child_edge,
        last_child_edge,
        tipname_to_id,
        ordered_tip_labels,
        node_labels,
    )
end

function _multi2di_node_data(tree::CompactTree)
    parent_of_node = copy(tree.parent_of_node)
    children = [copy(ch) for ch in tree.children]
    edge_length_to_node = zeros(Float64, tree.nnodes)
    @inbounds for edge in 1:tree.nedges
        edge_length_to_node[Int(tree.child_of_edge[edge])] = tree.edge_length[edge]
    end
    is_tip = copy(tree.is_tip)
    tip_labels = copy(tree.tip_labels)
    node_labels = copy(tree.node_labels)
    return parent_of_node, children, edge_length_to_node, is_tip, tip_labels, node_labels
end

function _multi2di_children_items(children_ids::Vector{Int32}, edge_length_to_node::Vector{Float64})
    items = Vector{_Multi2DiItem}(undef, length(children_ids))
    @inbounds for i in eachindex(children_ids)
        child = children_ids[i]
        items[i] = _Multi2DiItem(root = child, length = edge_length_to_node[Int(child)])
    end
    return items
end

function _multi2di_new_internal!(parent_of_node::Vector{Int32}, children::Vector{Vector{Int32}}, edge_length_to_node::Vector{Float64}, is_tip::BitVector, node_labels::Vector{String})
    push!(parent_of_node, Int32(0))
    push!(children, Int32[])
    push!(edge_length_to_node, 0.0)
    push!(is_tip, false)
    push!(node_labels, "")
    return Int32(length(parent_of_node))
end

function _multi2di_attach_child!(parent::Int32, child::Int32, length::Float64, parent_of_node::Vector{Int32}, children::Vector{Vector{Int32}}, edge_length_to_node::Vector{Float64})
    parent_of_node[Int(child)] = parent
    edge_length_to_node[Int(child)] = length
    push!(children[Int(parent)], child)
    return nothing
end

function _multi2di_resolve_ladder!(node::Int32, items::Vector{_Multi2DiItem}, parent_of_node::Vector{Int32}, children::Vector{Vector{Int32}}, edge_length_to_node::Vector{Float64}, is_tip::BitVector, node_labels::Vector{String})
    children[Int(node)] = Int32[]
    current = items[end]
    for idx in (length(items) - 1):-1:2
        newnode = _multi2di_new_internal!(parent_of_node, children, edge_length_to_node, is_tip, node_labels)
        _multi2di_attach_child!(newnode, items[idx].root, items[idx].length, parent_of_node, children, edge_length_to_node)
        _multi2di_attach_child!(newnode, current.root, current.length, parent_of_node, children, edge_length_to_node)
        current = _Multi2DiItem(root = newnode, length = 0.0)
    end
    _multi2di_attach_child!(node, items[1].root, items[1].length, parent_of_node, children, edge_length_to_node)
    _multi2di_attach_child!(node, current.root, current.length, parent_of_node, children, edge_length_to_node)
    return nothing
end

function _multi2di_resolve_random_pairs!(node::Int32, items::Vector{_Multi2DiItem}, parent_of_node::Vector{Int32}, children::Vector{Vector{Int32}}, edge_length_to_node::Vector{Float64}, is_tip::BitVector, node_labels::Vector{String}, rng::AbstractRNG)
    pool = copy(items)
    while length(pool) > 2
        i = rand(rng, 1:length(pool))
        j = rand(rng, 1:(length(pool) - 1))
        j >= i && (j += 1)
        a, b = minmax(i, j)
        item_a = pool[a]
        item_b = pool[b]
        deleteat!(pool, b)
        deleteat!(pool, a)
        newnode = _multi2di_new_internal!(parent_of_node, children, edge_length_to_node, is_tip, node_labels)
        _multi2di_attach_child!(newnode, item_a.root, item_a.length, parent_of_node, children, edge_length_to_node)
        _multi2di_attach_child!(newnode, item_b.root, item_b.length, parent_of_node, children, edge_length_to_node)
        push!(pool, _Multi2DiItem(root = newnode, length = 0.0))
    end
    children[Int(node)] = Int32[]
    _multi2di_attach_child!(node, pool[1].root, pool[1].length, parent_of_node, children, edge_length_to_node)
    _multi2di_attach_child!(node, pool[2].root, pool[2].length, parent_of_node, children, edge_length_to_node)
    return nothing
end

function multi2di(tree::CompactTree; random::Bool = true, equiprob::Bool = true, rng::AbstractRNG = Random.default_rng())
    parent_of_node, children, edge_length_to_node, is_tip, tip_labels, node_labels = _multi2di_node_data(tree)
    targets = [Int32(node) for node in tree.preorder if !tree.is_tip[Int(node)] && length(children[Int(node)]) > 2]
    isempty(targets) && return tree

    for node in targets
        items = _multi2di_children_items(children[Int(node)], edge_length_to_node)
        if random
            if equiprob
                _multi2di_resolve_random_pairs!(node, items, parent_of_node, children, edge_length_to_node, is_tip, node_labels, rng)
            else
                shuffle!(rng, items)
                _multi2di_resolve_ladder!(node, items, parent_of_node, children, edge_length_to_node, is_tip, node_labels)
            end
        else
            _multi2di_resolve_ladder!(node, items, parent_of_node, children, edge_length_to_node, is_tip, node_labels)
        end
    end

    return _compact_tree_from_topology(parent_of_node, children, edge_length_to_node, is_tip, tip_labels, node_labels)
end

function multi2di(tree::T; kwargs...) where {T <: NewickTree.Node}
    return from_compact_tree(multi2di(to_compact_tree(tree); kwargs...))
end
