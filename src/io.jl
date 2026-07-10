function load_newick_tree(path::AbstractString)
    isfile(path) || throw(ArgumentError("Tree file does not exist: $path"))
    newick = String(strip(read(path, String)))
    endswith(newick, ';') || throw(ArgumentError("Malformed Newick file: missing trailing semicolon"))
    try
        return NewickTree.readnw(newick)
    catch first_error
        cleaned = String(replace(newick, r"\)([^():;,]+):" => "):"))
        cleaned == newick && rethrow(first_error)
        try
            return NewickTree.readnw(cleaned)
        catch
            rethrow(first_error)
        end
    end
end

@inline function _safe_node_children(node::T) where {T}
    if !isdefined(node, :children)
        return T[]
    end
    raw_children = getfield(node, :children)
    child_nodes = T[]
    for i in eachindex(raw_children)
        if isassigned(raw_children, i)
            push!(child_nodes, raw_children[i])
        end
    end
    return child_nodes
end

function to_compact_tree(tree::T) where {T <: NewickTree.Node}
    root_node = tree
    while !NewickTree.isroot(root_node)
        root_node = root_node.parent
    end

    parent_of_node = Int32[]
    children = Vector{Int32}[]
    node_edge_length = Float64[]
    is_tip = BitVector()
    node_labels = String[]
    node_to_idx = IdDict{T, Int32}()

    stack = Tuple{T, Int32}[(root_node, 0)]
    while !isempty(stack)
        node, parent_idx = pop!(stack)
        haskey(node_to_idx, node) && continue

        idx = Int32(length(parent_of_node) + 1)
        node_to_idx[node] = idx
        push!(parent_of_node, parent_idx)
        push!(children, Int32[])
        push!(node_edge_length, NewickTree.isroot(node) ? 0.0 : Float64(NewickTree.distance(node)))
        push!(node_labels, String(NewickTree.name(node)))

        child_nodes = _safe_node_children(node)
        push!(is_tip, isempty(child_nodes))

        if parent_idx != 0
            push!(children[parent_idx], idx)
        end

        for child in reverse(child_nodes)
            push!(stack, (child, idx))
        end
    end

    nnodes = length(parent_of_node)
    root = get(node_to_idx, root_node, Int32(0))
    root == 0 && throw(ErrorException("Failed to identify root node while building CompactTree"))

    dist_from_root = zeros(Float64, nnodes)
    preorder = Int32[]
    stack2 = Int32[root]
    while !isempty(stack2)
        node = pop!(stack2)
        push!(preorder, node)
        for child in Iterators.reverse(children[node])
            dist_from_root[child] = dist_from_root[node] + node_edge_length[child]
            push!(stack2, child)
        end
    end

    postorder = Int32[]
    postorder_internal = Int32[]
    function dfs(node::Int32)
        for child in children[node]
            dfs(child)
        end
        push!(postorder, node)
        !is_tip[node] && push!(postorder_internal, node)
        return nothing
    end
    dfs(root)

    tip_ids = Int32[findall(is_tip)...]
    ntips = length(tip_ids)
    nedges = nnodes - 1
    tip_labels = [node_labels[idx] for idx in tip_ids]
    tipname_to_id = Dict{String, Int32}(tip_labels[i] => tip_ids[i] for i in eachindex(tip_ids))

    parent_of_edge = Vector{Int32}(undef, nedges)
    child_of_edge = Vector{Int32}(undef, nedges)
    edge_length = Vector{Float64}(undef, nedges)
    first_child_edge = fill(Int32(0), nnodes)
    last_child_edge = fill(Int32(0), nnodes)

    edge_idx = Int32(1)
    for parent in Int32.(1:nnodes)
        for child in children[parent]
            parent_of_edge[edge_idx] = parent
            child_of_edge[edge_idx] = child
            edge_length[edge_idx] = node_edge_length[child]
            if first_child_edge[parent] == 0
                first_child_edge[parent] = edge_idx
            end
            last_child_edge[parent] = edge_idx
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
        tip_labels,
        node_labels,
    )
end

preprocess_tree(tree::T) where {T <: NewickTree.Node} = to_compact_tree(tree)

function _build_children_with_lengths(tree::CompactTree)
    children_dl = Dict{Int, Vector{Tuple{Int, Float64}}}()
    for e in 1:tree.nedges
        p = Int(tree.parent_of_edge[e])
        c = Int(tree.child_of_edge[e])
        d = tree.edge_length[e]
        list = get!(children_dl, p) do
            Tuple{Int, Float64}[]
        end
        push!(list, (c, d))
    end
    return children_dl
end

function _compact_tip_label(tree::CompactTree, node::Int)
    idx = findfirst(==(Int32(node)), tree.tip_ids)
    idx === nothing && return ""
    return tree.tip_labels[idx]
end

function _compact_internal_label(tree::CompactTree, node::Int)
    label = tree.node_labels[node]
    return isempty(label) ? "" : label
end

function _build_nw_node_with_distance(tree::CompactTree, node::Int, d::Float64, children_dl::Dict{Int, Vector{Tuple{Int, Float64}}})
    has_children = haskey(children_dl, node)
    label = tree.is_tip[node] ? _compact_tip_label(tree, node) : _compact_internal_label(tree, node)
    if !has_children
        return NewickTree.Node(UInt32(node), NewickData(d = d, n = label))
    end
    nw = NewickTree.Node(UInt32(node), NewickData(d = d, n = label))
    for (child, cd) in children_dl[node]
        push!(nw, _build_nw_node_with_distance(tree, child, cd, children_dl))
    end
    return nw
end

function _build_nw_node(tree::CompactTree, node::Int, children_dl::Dict{Int, Vector{Tuple{Int, Float64}}})
    has_children = haskey(children_dl, node)
    label = tree.is_tip[node] ? _compact_tip_label(tree, node) : _compact_internal_label(tree, node)
    if !has_children
        return NewickTree.Node(UInt32(node), NewickData(d = 0.0, n = label))
    end
    nw = NewickTree.Node(UInt32(node), NewickData(d = 0.0, n = label))
    for (child, d) in children_dl[node]
        push!(nw, _build_nw_node_with_distance(tree, child, d, children_dl))
    end
    return nw
end

function from_compact_tree(tree::CompactTree)
    children_dl = _build_children_with_lengths(tree)
    return _build_nw_node(tree, Int(tree.root), children_dl)
end

function save_newick_tree(path::AbstractString, tree::CompactTree; append::Bool = false)
    io = IOBuffer()
    NewickTree.writenw(io, from_compact_tree(tree))
    return _write_newick_file(path, String(take!(io)), append)
end

function save_newick_tree(path::AbstractString, nw::NewickTree.Node; append::Bool = false)
    io = IOBuffer()
    NewickTree.writenw(io, nw)
    return _write_newick_file(path, String(take!(io)), append)
end

function save_newick_tree(path::AbstractString, newick::AbstractString; append::Bool = false)
    content = String(newick)
    endswith(strip(content), ';') || (content = string(strip(content), ';'))
    return _write_newick_file(path, content, append)
end

function _write_newick_file(path::AbstractString, content::AbstractString, append::Bool)
    open(path, append ? "a" : "w") do io
        write(io, content)
        write(io, '\n')
    end
    return String(content)
end
