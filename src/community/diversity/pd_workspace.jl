# Presence counts on a site tree describe every original edge on a contracted
# path. Rooted variants also include the original root-to-LCA path; its length
# is summed by ordered path segments, never by subtracting two root heights.
mutable struct _PDWorkspace{G}
    geometry::G
    index::Union{Nothing, _SpeciesPathIndex}
    site::Union{Nothing, _SpeciesSiteWorkspace}
    counts::Vector{Int}
    nodes::Vector{Int32}
end

function _PDWorkspace(tree)
    lengths = _comdist_lengths(tree)
    geometry = (; factors = ones(tree.nnodes), offsets = lengths,
        down_offsets = lengths, heights = Float64[], scaled = false)
    return _PDWorkspace(geometry, nothing, nothing, Int[], Int32[])
end

@inline function _pd_edge_selected(count, richness, variant)
    variant === :tbl && return count > 0
    variant === :sbl && return count > 1
    variant === :ueh && return count == 1
    return 0 < count < richness # MST; TIP uses original terminal edges.
end

function _pd_site!(workspace, tree, tips, variant)
    r = length(tips)
    r == 0 && return 0.0
    if variant === :tip
        value = 0.0
        for tip in tips
            value += workspace.geometry.offsets[tree.tip_ids[tip]]
        end
        return value
    end
    # Sequential counts win for densely occupied sites, without virtual-tree
    # sorting or indexing. All arrays are retained across bootstrap draws.
    if !_species_site_virtual(tree.nnodes, r) && r > 1
        counts = workspace.counts
        resize!(counts, tree.nnodes)
        fill!(counts, 0)
        for tip in tips
            counts[tree.tip_ids[tip]] = 1
        end
        for node in tree.postorder
            parent = tree.parent_of_node[node]
            parent != 0 && (counts[parent] += counts[node])
        end
        value = 0.0
        for edge in 1:tree.nedges
            _pd_edge_selected(counts[tree.child_of_edge[edge]], r, variant) &&
                (value += tree.edge_length[edge])
        end
        return value
    end
    if workspace.index === nothing
        workspace.index = _SpeciesPathIndex(tree, workspace.geometry)
        workspace.site = _SpeciesSiteWorkspace(tree.nnodes)
    end
    nodes = workspace.nodes
    resize!(nodes, r)
    for (i, tip) in enumerate(tips)
        nodes[i] = tree.tip_ids[tip]
    end
    plan = (; tree, geometry = workspace.geometry, index = workspace.index, strategy = :virtual)
    topology, geometry, local_tips = _species_site!(workspace.site, plan, nodes)
    counts = workspace.counts
    resize!(counts, topology.nnodes)
    fill!(counts, 0)
    for tip in local_tips
        counts[tip] = 1
    end
    value = 0.0
    for node in topology.postorder
        parent = topology.parent_of_node[node]
        parent == 0 && continue
        _pd_edge_selected(counts[node], r, variant) && (value += geometry.offsets[node])
        counts[parent] += counts[node]
    end
    if _pd_edge_selected(r, r, variant)
        ancestor = workspace.site.nodes[1]
        value += _species_path(workspace.index, tree, tree.root, ancestor)[2]
    end
    return value
end

function _pd_tip_positions(tree, tips)
    values = tips isa AbstractString ? [tips] : collect(tips)
    # Resolve once; aliases such as position 1 and its name are duplicates too.
    positions = Int32[]
    inverse = nothing
    for value in values
        if value isa Integer && !(value isa Bool)
            1 <= value <= tree.ntips || throw(ArgumentError("tip position is outside the tree"))
            push!(positions, Int32(value))
        elseif value isa AbstractString
            node = get(tree.tipname_to_id, String(value), Int32(0))
            node != 0 && tree.is_tip[node] || throw(ArgumentError("unknown tree tip: $value"))
            if inverse === nothing
                inverse = zeros(Int32, tree.nnodes)
                for (i, tip) in enumerate(tree.tip_ids)
                    inverse[tip] = i
                end
            end
            push!(positions, inverse[node])
        else
            throw(ArgumentError("tips must contain integer positions or tip names"))
        end
    end
    sort!(positions)
    allunique(positions) || throw(ArgumentError("tips must not contain duplicates"))
    return positions
end

function _pd_integer(value, name, minimum)
    !(value isa Bool) && minimum <= value <= typemax(Int) ||
        throw(ArgumentError("$name must be an integer between $minimum and typemax(Int)"))
    return Int(value)
end

function _pd_probabilities(tree, tip_weights)
    weights = tip_weights === nothing ? ones(tree.ntips) : Float64.(tip_weights)
    length(weights) == tree.ntips || throw(ArgumentError("tip_weights must match tree tips"))
    all(isfinite, weights) && all(>=(0), weights) && any(>(0), weights) ||
        throw(ArgumentError("tip_weights must be finite and non-negative with positive sum"))
    total = sum(weights)
    if !isfinite(total)
        # Keep normal-range keys and RNG consumption unchanged. Scaling only
        # repairs overflowing sums, which previously made every probability 0.
        weights ./= maximum(weights)
        total = sum(weights)
    end
    return weights ./ total
end

# Exponential-key sampling without replacement, with the original random
# draws and index-stable ties. Keep min(r,p-r) heap entries: the best r keys
# for a sparse draw, or the worst p-r for a dense draw. No full sort/key vector.
struct _PDSamplingWorkspace
    keys::Vector{Float64}
    heap::Vector{Int32}
    selected::Vector{Int32}
    excluded::BitVector
end
_PDSamplingWorkspace() = _PDSamplingWorkspace(Float64[], Int32[], Int32[], falses(0))

@inline _pd_key_better(a, i, b, j) = a > b || (a == b && i < j)
@inline _pd_heap_before(a, i, b, j, best) = best ? _pd_key_better(b, j, a, i) : _pd_key_better(a, i, b, j)

function _pd_sample!(workspace, probabilities, r, rng)
    p = length(probabilities)
    best = r <= p ÷ 2
    k = min(r, p - r)
    keys, heap, selected = workspace.keys, workspace.heap, workspace.selected
    resize!(keys, k)
    resize!(heap, k)
    resize!(selected, r)
    used = 0
    for tip in 1:p
        probability = probabilities[tip]
        key = probability > 0 ? log(rand(rng)) / probability : -Inf
        k == 0 && continue # Still consume precisely the original RNG draws.
        if used < k
            used += 1
            position = used
            while position > 1
                parent = position ÷ 2
                _pd_heap_before(key, tip, keys[parent], heap[parent], best) || break
                keys[position], heap[position] = keys[parent], heap[parent]
                position = parent
            end
            keys[position], heap[position] = key, tip
        elseif best ? _pd_key_better(key, tip, keys[1], heap[1]) : _pd_key_better(keys[1], heap[1], key, tip)
            position = 1
            while 2position <= k
                child = 2position
                if child < k && _pd_heap_before(keys[child + 1], heap[child + 1], keys[child], heap[child], best)
                    child += 1
                end
                _pd_heap_before(keys[child], heap[child], key, tip, best) || break
                keys[position], heap[position] = keys[child], heap[child]
                position = child
            end
            keys[position], heap[position] = key, tip
        end
    end
    if best
        copyto!(selected, heap)
        sort!(selected)
    else
        resize!(workspace.excluded, p)
        fill!(workspace.excluded, false)
        for tip in heap
            workspace.excluded[tip] = true
        end
        position = 0
        for tip in 1:p
            if !workspace.excluded[tip]
                selected[position += 1] = tip
            end
        end
    end
    return selected
end
