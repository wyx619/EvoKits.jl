# Contract only internal nodes whose incoming edge is exactly zero. The root
# and every tip survive; preorder attachment preserves the order of children,
# and increasing old IDs preserve the existing tip/column order.
function _timeslice_collapse(tree::CompactTree, lengths::Vector{Float64})
    incoming = zeros(Float64, tree.nnodes)
    for edge in eachindex(lengths)
        incoming[Int(tree.child_of_edge[edge])] = lengths[edge]
    end
    retained = tree.is_tip .| (incoming .> 0)
    retained[Int(tree.root)] = true
    node_map = zeros(Int32, tree.nnodes)
    count = 0
    for node in 1:tree.nnodes
        if retained[node]
            count += 1
            node_map[node] = Int32(count)
        end
    end
    parents = zeros(Int32, count)
    children = [Int32[] for _ in 1:count]
    new_lengths = zeros(Float64, count)
    is_tip = falses(count)
    labels = Vector{String}(undef, count)
    owner = zeros(Int32, tree.nnodes)
    for node32 in tree.preorder
        node = Int(node32)
        old_parent = Int(tree.parent_of_node[node])
        parent = old_parent == 0 ? Int32(0) : owner[old_parent]
        if retained[node]
            new_node = node_map[node]
            owner[node] = new_node
            parents[new_node] = parent
            parent != 0 && push!(children[parent], new_node)
            new_lengths[new_node] = incoming[node]
            is_tip[new_node] = tree.is_tip[node]
            labels[new_node] = tree.node_labels[node]
        else
            owner[node] = parent
        end
    end
    return _compact_tree_from_topology(parents, children, new_lengths,
        is_tip, copy(tree.tip_labels), labels)
end

"""
    timeslice(tree::CompactTree, depth; collapse=false, atol=0.0, rtol=1e-8)

Remove the most recent `depth` units of branch history from a dated ultrametric
tree. Depth is measured backwards from the present in the tree's branch-length
units and must lie in [0, tree height]. Return an independent `CompactTree`;
retain every species and the original tip/column order, including zero-length
terminal edges. No input is modified or automatically dated.

An edge starting at parent age `a` retains `min(length, max(a-depth, 0))`.
Branches entirely older than the cut are unchanged, crossing edges are
shortened, and branches starting exactly at or younger than the cut become
zero. Depth zero preserves lengths; depth equal to the height makes all edges
zero. Ages use longest descendant paths, avoiding H-root_depth cancellation
below a long shared stem. Validation and atol/rtol match [`get_clades`](@ref);
tolerances do not shift the cut or decide which edges are zero.

By default topology, node IDs, edge order and node labels are retained, and
root distances are recomputed. With `collapse=true`, contract exactly zero
internal edges (including pre-existing ones), retain the root and all tips,
and rebuild IDs and traversal caches; labels on contracted nodes disappear.
Positive edges, however short, are never contracted.

The result can be used with the original species-level `CommunityData` for
`pd`, `phylo_endemism`, `phylobeta` and regionalization. Prepare a new
`CommunityContext` for each sliced tree; contexts are snapshots. This changes
branch lengths rather than merging community columns as `collapse_range`
does. Rooted PD/PE and length-based beta are invariant to zero-edge collapse;
edge-indexed diagnostics and unweighted topology statistics need not be.

Time and additional storage are O(nodes+tips), with no pairwise tip matrix.
For several absolute time depths, slice the original tree separately; slicing
an already sliced tree removes additional history from its new present.
"""
function timeslice(tree::CompactTree, depth; collapse::Bool = false,
    atol = 0.0, rtol = 1e-8,
)
    selected = _clade_time_value(depth, "depth")
    abs_tolerance = _clade_time_value(atol, "atol")
    rel_tolerance = _clade_time_value(rtol, "rtol")
    length(tree.tip_labels) == tree.ntips && allunique(tree.tip_labels) &&
        all(label -> !isempty(strip(label)), tree.tip_labels) ||
        throw(ArgumentError("tree tip labels must be unique and nonempty"))
    ages = _clade_node_ages(tree, abs_tolerance, rel_tolerance)
    selected <= ages[Int(tree.root)] ||
        throw(ArgumentError("depth must not exceed the tree height"))
    lengths = [min(tree.edge_length[edge],
        max(ages[Int(tree.parent_of_edge[edge])] - selected, 0.0))
        for edge in eachindex(tree.edge_length)]
    collapse && return _timeslice_collapse(tree, lengths)

    depths = zeros(Float64, tree.nnodes)
    for node32 in tree.preorder
        node = Int(node32)
        tree.is_tip[node] && continue
        for edge in Int(tree.first_child_edge[node]):Int(tree.last_child_edge[node])
            child = Int(tree.child_of_edge[edge])
            depths[child] = depths[node] + lengths[edge]
            isfinite(depths[child]) || throw(ArgumentError("root-to-tip lengths must remain finite"))
        end
    end
    return _compact_tree_with_lengths(tree, lengths, depths)
end
