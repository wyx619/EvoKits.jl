# Build the same minimal connecting tree used by LCA-based virtual trees,
# using one linear postorder pass instead of constructing a binary LCA index.
# Only observed tips and nodes joining >=2 active child subtrees survive.
# Original edges are summed along contracted paths, never by root subtraction.
function _comdist_virtual(rows, tree)
    representative = zeros(Int32, tree.nnodes)
    for nodes in rows.nodes, node in nodes
        representative[node] = -1
    end
    parents, lengths, children = Int32[], Float64[], Vector{Int32}[]
    for node in tree.postorder
        if representative[node] == -1
            push!(parents, 0)
            push!(lengths, 0.0)
            push!(children, Int32[])
            representative[node] = length(parents)
            continue
        end
        first = tree.first_child_edge[node]
        first == 0 && continue
        count, only_child = 0, Int32(0)
        for edge in first:tree.last_child_edge[node]
            child = representative[tree.child_of_edge[edge]]
            child == 0 && continue
            lengths[child] += tree.edge_length[edge]
            count += 1
            only_child = child
        end
        if count == 1
            representative[node] = only_child
        elseif count > 1
            push!(parents, 0)
            push!(lengths, 0.0)
            push!(children, Int32[])
            parent = Int32(length(parents))
            representative[node] = parent
            for edge in first:tree.last_child_edge[node]
                child = representative[tree.child_of_edge[edge]]
                child == 0 && continue
                parents[child] = parent
                push!(children[parent], child)
            end
        end
    end
    m = length(parents)
    # Discard the shared stem above the observed tips' LCA.
    m > 0 && (lengths[m] = 0.0)
    nodes = [Int32[representative[node] for node in row] for row in rows.nodes]
    topology = (; nnodes = m, root = Int32(m), parent_of_node = parents,
        children, postorder = 1:m, preorder = m:-1:1)
    return merge(rows, (; nodes, lengths)), topology
end

# Private strategy override supports numerical cross-checks and benchmarks;
# the public API automatically selects a route and keeps its return contract.
function _comdist_use_virtual(rows, tree, targets)
    tree.nnodes >= 256 && targets >= 1 || return false
    # Stop as soon as occupancy rules out compression; dense repeated records
    # need not all be scanned merely to select the full-tree route.
    budget = tree.nnodes ÷ (targets == 1 ? 16 : 4)
    max_tips = (budget + 1) ÷ 2
    seen = falses(tree.nnodes)
    u = 0
    for nodes in rows.nodes, node in nodes
        if !seen[node]
            seen[node] = true
            u += 1
            u > max_tips && return false
        end
    end
    # The connecting tree has at most 2u-1 nodes for u observed tips.
    return true
end

function _comdist_selected(rows, tree, selected, order)
    result = fill(NaN, length(selected))
    isempty(order) && return result
    workspace = _ComdistWorkspace(tree.nnodes)
    last_target, last_source, value = 0, 0, NaN
    for output in order
        target, source = minmax(selected[output]...)
        if isempty(rows.nodes[target]) || isempty(rows.nodes[source])
            continue
        end
        if target != last_target
            _comdist_target!(workspace, tree, rows.nodes[target], rows.weights[target], rows.lengths)
            last_target, last_source = target, 0
        end
        if source != last_source
            value = 0.0
            for (node, weight) in zip(rows.nodes[source], rows.weights[source])
                value += weight * (workspace.down[node] + workspace.outside_distance[node])
            end
            last_source = source
        end
        result[output] = value
    end
    return result
end

function _comdist_adaptive(rows, tree, selected; strategy = :auto)
    strategy in (:auto, :full, :virtual) || throw(ArgumentError("unknown private comdist strategy"))
    if selected === nothing
        targets = count(!isempty, rows.nodes)
        targets == 0 && return fill(NaN, length(rows.nodes), length(rows.nodes))
        order = Int[]
    else
        isempty(selected) && return Float64[]
        order = sortperm(selected; by = pair -> minmax(pair...))
        targets, last_target = 0, 0
        used = falses(length(rows.nodes))
        for output in order
            target, source = minmax(selected[output]...)
            (isempty(rows.nodes[target]) || isempty(rows.nodes[source])) && continue
            used[target] = used[source] = true
            if target != last_target
                targets += 1
                last_target = target
            end
        end
        targets == 0 && return fill(NaN, length(selected))
        # Endpoints paired only with empty sites cannot affect any value.
        # Exclude their tips from the connecting tree without mutating rows.
        if any(i -> !used[i] && !isempty(rows.nodes[i]), eachindex(used))
            rows = merge(rows, (; nodes = [used[i] ? row : Int32[] for (i, row) in enumerate(rows.nodes)]))
        end
    end
    virtual = strategy === :virtual || (strategy === :auto && _comdist_use_virtual(rows, tree, targets))
    if virtual
        rows, tree = _comdist_virtual(rows, tree)
    elseif isempty(rows.lengths)
        rows = merge(rows, (; lengths = _comdist_lengths(tree)))
    end
    return selected === nothing ? _comdist(rows, tree) : _comdist_selected(rows, tree, selected, order)
end
