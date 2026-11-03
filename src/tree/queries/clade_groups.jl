"""
    CladeGroupingResult

Time-slice partition of a complete tree returned by [`get_clades`](@ref).
`tip_labels` and `membership` follow input tip order; `group_labels`,
`group_nodes` (CompactTree node IDs), and `group_sizes` follow first-tip order.
`cut` is age before the present, `boundary` is `:older` or `:younger`, and
`tree_height` is the maximum root-to-tip length. `requested_k` is nothing for
an explicit cut; `actual_k` records the achieved number. Arrays are owned.
"""
struct CladeGroupingResult
    tip_labels::Vector{String}
    membership::Vector{Int}
    group_labels::Vector{String}
    group_nodes::Vector{Int}
    group_sizes::Vector{Int}
    cut::Float64
    boundary::Symbol
    tree_height::Float64
    requested_k::Union{Nothing, Int}
    actual_k::Int
end

function _clade_time_value(value, name)
    value isa Real && !(value isa Bool) && isfinite(value) && value >= 0 ||
        throw(ArgumentError("$name must be a finite non-negative real number"))
    converted = Float64(value)
    isfinite(converted) && (value == 0 || converted > 0) ||
        throw(ArgumentError("$name must be representable as finite Float64 without underflow to zero"))
    return converted
end

# Ages are longest remaining paths, accumulated from tips. Unlike H-depth,
# this retains short young branches below a long shared stem. The shortest
# remaining path at the root checks ultrametricity without cached depths.
function _clade_node_ages(tree, atol, rtol)
    tree.ntips > 0 || throw(ArgumentError("tree must contain at least one tip"))
    ages, shortest = zeros(tree.nnodes), zeros(tree.nnodes)
    for edge in eachindex(tree.edge_length)
        value = tree.edge_length[edge]
        isfinite(value) && value >= 0 || throw(ArgumentError("tree branch lengths must be finite and non-negative"))
    end
    for node32 in tree.postorder_internal
        node = Int(node32)
        lo, hi = Inf, 0.0
        for edge in Int(tree.first_child_edge[node]):Int(tree.last_child_edge[node])
            child = Int(tree.child_of_edge[edge])
            lo = min(lo, shortest[child] + tree.edge_length[edge])
            hi = max(hi, ages[child] + tree.edge_length[edge])
        end
        isfinite(lo) && isfinite(hi) || throw(ArgumentError("root-to-tip lengths must remain finite"))
        shortest[node], ages[node] = lo, hi
    end
    root = Int(tree.root)
    isapprox(shortest[root], ages[root]; atol, rtol) ||
        throw(ArgumentError("time-slice grouping requires an ultrametric tree; date or transform the tree explicitly"))
    return ages
end

function _clade_target_cut(tree, ages, k)
    k isa Integer && !(k isa Bool) && 1 <= k <= tree.ntips ||
        throw(ArgumentError("k must be an integer between 1 and the number of tips"))
    k = Int(k)
    k == 1 && return ages[Int(tree.root)], :older, k
    events = [Int(node) for node in tree.postorder_internal if length(tree.children[node]) > 1]
    sort!(events; by = node -> ages[node], rev = true)
    count, index = 1, 1
    while index <= length(events)
        age = ages[events[index]]
        while index <= length(events) && ages[events[index]] == age
            count += length(tree.children[events[index]]) - 1
            index += 1
        end
        count >= k && return age, :younger, k
    end
    throw(ArgumentError("tree does not contain enough branching events for k"))
end

"""
    get_clades(tree::CompactTree; cut=nothing, k=nothing, boundary=nothing,
               atol=0.0, rtol=1e-8)

Partition a complete dated ultrametric tree into disjoint lineage groups.
Supply exactly one of `cut` (age before the present, 0 <= cut <= tree height)
or target `k` (integer 1:ntips). Return a `CladeGroupingResult` with each tip's
group, source node IDs, labels, sizes and achieved cut metadata. Groups are
numbered by their first input tip; all tree tips are retained.

For explicit cuts, `boundary=:older` (default) keeps a node exactly at the
cut as one group; `:younger` splits it. All nodes at the same represented age
are treated together, including zero-length edges. At cut=0, older may retain
zero-age clades while younger gives singleton tips; at tree height, older
gives one group. No epsilon or tolerance is applied to boundary comparisons.

For target k, choose the oldest complete branching event reaching at least k
groups, splitting every simultaneous node together. Multifurcations or equal
ages can yield `actual_k > requested_k`; no tied event is split arbitrarily.
k=1 keeps the root. The result records the exact event age and boundary, so
`get_clades(tree; cut=result.cut, boundary=result.boundary)` reproduces it.
Do not specify boundary with k. Unary nodes do not increase the group count.

Edges must be finite and non-negative, and root-to-tip lengths must be finite.
Ultrametricity is checked by shortest/longest root-to-tip paths using atol/rtol;
accepted near-ultrametric trees use longest descendant path as each node age.
Tolerances only validate the tree, never merge event ages. No tree is changed,
pruned or automatically dated; arrays are independent. Explicit cuts take
O(nodes+tips) time/storage; target k additionally sorts branching events.
"""
function get_clades(tree::CompactTree; cut = nothing, k = nothing,
    boundary = nothing, atol = 0.0, rtol = 1e-8,
)
    (cut === nothing) != (k === nothing) || throw(ArgumentError("supply exactly one of cut or k"))
    abs_tolerance = _clade_time_value(atol, "atol")
    rel_tolerance = _clade_time_value(rtol, "rtol")
    labels = copy(tree.tip_labels)
    length(labels) == tree.ntips && allunique(labels) && all(label -> !isempty(strip(label)), labels) ||
        throw(ArgumentError("tree tip labels must be unique and nonempty"))
    ages = _clade_node_ages(tree, abs_tolerance, rel_tolerance)
    height = ages[Int(tree.root)]
    if k === nothing
        selected = _clade_time_value(cut, "cut")
        selected <= height || throw(ArgumentError("cut must not exceed the tree height"))
        side = boundary === nothing ? :older : boundary
        side in (:older, :younger) || throw(ArgumentError("boundary must be :older or :younger"))
        requested = nothing
    else
        boundary === nothing || throw(ArgumentError("boundary is determined by the complete branching event when k is supplied"))
        selected, side, requested = _clade_target_cut(tree, ages, k)
    end
    owners = zeros(Int, tree.nnodes)
    for node32 in tree.preorder
        node = Int(node32)
        parent = Int(tree.parent_of_node[node])
        if parent != 0 && owners[parent] != 0
            owners[node] = owners[parent]
        elseif tree.is_tip[node] || ages[node] < selected || (side === :older && ages[node] == selected)
            owners[node] = node
        end
    end
    membership, nodes, sizes = Vector{Int}(undef, tree.ntips), Int[], Int[]
    groups = zeros(Int, tree.nnodes)
    for (tip, node32) in enumerate(tree.tip_ids)
        owner = owners[Int(node32)]
        if groups[owner] == 0
            push!(nodes, owner)
            push!(sizes, 0)
            groups[owner] = length(nodes)
        end
        group = groups[owner]
        membership[tip] = group
        sizes[group] += 1
    end
    group_labels = ["clade_$i" for i in eachindex(nodes)]
    return CladeGroupingResult(labels, membership, group_labels, nodes, sizes,
        selected, side, height, requested, length(nodes))
end
