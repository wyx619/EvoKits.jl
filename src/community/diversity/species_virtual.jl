# An ordered path stores its factor and both affine offsets. For edges A
# above B, lifting applies A(B(x)); descending applies B(A(x)). Keep the
# offsets directly: 1-product(factors) would erase tiny positive contrasts.
@inline _species_join(a, b) = (a[1] * b[1],
    a[2] + a[1] * b[2], b[3] + b[1] * a[3])
const _SPECIES_IDENTITY = (1.0, 0.0, 0.0)

# One postorder pass contracts unobserved/unary paths. Original heights and
# edge transforms encode the shared root stem even when the working root is
# the observed LCA. Raw PSC uses these original heights and can be negative.
function _species_compress(nodes, tree, geometry)
    representative = zeros(Int32, tree.nnodes)
    for row in nodes, node in row
        representative[node] = -1
    end
    parents, children = Int32[], Vector{Int32}[]
    factors, offsets, down_offsets, heights = Float64[], Float64[], Float64[], Float64[]
    function append_node(node)
        push!(parents, 0)
        push!(children, Int32[])
        push!(factors, 1.0)
        push!(offsets, 0.0)
        push!(down_offsets, 0.0)
        !isempty(geometry.heights) && push!(heights, geometry.heights[node])
        return Int32(length(parents))
    end
    for node in tree.postorder
        if representative[node] == -1
            representative[node] = append_node(node)
            continue
        end
        count, only_child = 0, Int32(0)
        for original in tree.children[node]
            child = representative[original]
            child == 0 && continue
            f, up, down = _species_join(
                (geometry.factors[original], geometry.offsets[original], geometry.down_offsets[original]),
                (factors[child], offsets[child], down_offsets[child]))
            factors[child], offsets[child], down_offsets[child] = f, up, down
            count += 1
            only_child = child
        end
        if count == 1
            representative[node] = only_child
        elseif count > 1
            parent = append_node(node)
            representative[node] = parent
            for original in tree.children[node]
                child = representative[original]
                child == 0 && continue
                parents[child] = parent
                push!(children[parent], child)
            end
        end
    end
    m = length(parents)
    if m > 0
        factors[m], offsets[m], down_offsets[m] = _SPECIES_IDENTITY
    end
    mapped = [Int32[representative[node] for node in row] for row in nodes]
    topology = (; nnodes = m, root = Int32(m), parent_of_node = parents,
        children, postorder = 1:m, preorder = m:-1:1, tip_ids = Base.OneTo(m))
    return mapped, topology, (; factors, offsets, down_offsets, heights, scaled = geometry.scaled)
end

# Heavy paths plus an ordered segment tree use linear storage. They answer
# LCA in O(log v) and stable affine path queries in O(log^2 v), without an
# O(v log v) binary-lifting table or root-depth subtraction.
struct _SpeciesPathIndex
    head::Vector{Int32}
    position::Vector{Int32}
    sizes::Vector{Int32}
    depth::Vector{Int32}
    paths::Vector{NTuple{3, Float64}}
    base::Int
end

function _SpeciesPathIndex(tree, geometry)
    v = tree.nnodes
    sizes, heavy = ones(Int32, v), zeros(Int32, v)
    depth = zeros(Int32, v)
    for node in tree.postorder
        for child in tree.children[node]
            sizes[node] += sizes[child]
            (heavy[node] == 0 || sizes[child] > sizes[heavy[node]]) && (heavy[node] = child)
        end
    end
    for node in tree.preorder
        parent = tree.parent_of_node[node]
        parent != 0 && (depth[node] = depth[parent] + 1)
    end
    head, position = zeros(Int32, v), zeros(Int32, v)
    stack, rank = Int32[tree.root], Int32(0)
    while !isempty(stack)
        first = node = pop!(stack)
        while node != 0
            head[node], position[node] = first, (rank += 1)
            for child in Iterators.reverse(tree.children[node])
                child != heavy[node] && push!(stack, child)
            end
            node = heavy[node]
        end
    end
    base = nextpow(2, max(1, v))
    paths = fill(_SPECIES_IDENTITY, 2base)
    for node in tree.preorder
        paths[base + position[node] - 1] =
            (geometry.factors[node], geometry.offsets[node], geometry.down_offsets[node])
    end
    for i in (base - 1):-1:1
        paths[i] = _species_join(paths[2i], paths[2i + 1])
    end
    return _SpeciesPathIndex(head, position, sizes, depth, paths, base)
end

function _species_lca(index, tree, a, b)
    while index.head[a] != index.head[b]
        if index.depth[index.head[a]] > index.depth[index.head[b]]
            a = tree.parent_of_node[index.head[a]]
        else
            b = tree.parent_of_node[index.head[b]]
        end
    end
    return index.depth[a] < index.depth[b] ? a : b
end

@inline _species_ancestor(index, a, b) =
    index.position[a] <= index.position[b] < index.position[a] + index.sizes[a]

function _species_range(index, left, right)
    left += index.base - 1
    right += index.base - 1
    before = after = _SPECIES_IDENTITY
    while left <= right
        if isodd(left)
            before = _species_join(before, index.paths[left])
            left += 1
        end
        if iseven(right)
            after = _species_join(index.paths[right], after)
            right -= 1
        end
        left >>= 1
        right >>= 1
    end
    return _species_join(before, after)
end

function _species_path(index, tree, ancestor, node)
    result = _SPECIES_IDENTITY
    while index.head[ancestor] != index.head[node]
        head = index.head[node]
        result = _species_join(_species_range(index, index.position[head], index.position[node]), result)
        node = tree.parent_of_node[head]
    end
    return _species_join(_species_range(index, index.position[ancestor] + 1, index.position[node]), result)
end

_species_site_virtual(v, richness) = richness > 1 && v >= 256 && 2richness - 1 <= v ÷ 32

function _species_plan(rows, tree, geometry; strategy = :auto)
    strategy in (:auto, :full, :virtual) || throw(ArgumentError("unknown private species strategy"))
    targets = count(row -> length(row) > 1, rows)
    compress = strategy === :virtual || (strategy === :auto &&
        _comdist_use_virtual((; nodes = Iterators.filter(row -> length(row) > 1, rows)), tree, targets))
    if !compress && (strategy === :full ||
        !any(row -> _species_site_virtual(tree.nnodes, length(row)), rows))
        # Dense collections keep the original records/topology, with no node
        # remapping, copies or index. The established sequential kernel wins here.
        return (; rows, tree, geometry, index = nothing, strategy)
    end
    # Empty/singleton rows need no computation and must not enlarge the tree.
    nodes = [length(row) > 1 ? Int32[tree.tip_ids[tip] for tip in row] : Int32[] for row in rows]
    if compress
        nodes, topology, geometry = _species_compress(nodes, tree, geometry)
    else
        topology = (; nnodes = tree.nnodes, root = tree.root,
            parent_of_node = tree.parent_of_node, children = tree.children,
            postorder = tree.postorder, preorder = tree.preorder, tip_ids = Base.OneTo(tree.nnodes))
    end
    local_tree = strategy === :virtual ? any(!isempty, nodes) :
        strategy !== :full && any(row -> _species_site_virtual(topology.nnodes, length(row)), nodes)
    index = local_tree ? _SpeciesPathIndex(topology, geometry) : nothing
    return (; rows = nodes, tree = topology, geometry, index, strategy)
end

struct _SpeciesChildren
    first::Vector{Int32}
    next::Vector{Int32}
end
struct _SpeciesChildIterator
    node::Int32
    next::Vector{Int32}
end
Base.getindex(children::_SpeciesChildren, node) = _SpeciesChildIterator(children.first[node], children.next)
Base.iterate(children::_SpeciesChildIterator, node = children.node) =
    node == 0 ? nothing : (node, children.next[node])

# All site-local containers grow to the largest required working tree and are
# reused. Mapped tips retain abundance order even though virtual nodes are sorted.
struct _SpeciesSiteWorkspace
    selected::Vector{Int32}
    nodes::Vector{Int32}
    parents::Vector{Int32}
    stack::Vector{Int32}
    local_index::Vector{Int32}
    tips::Vector{Int32}
    children::_SpeciesChildren
    factors::Vector{Float64}
    offsets::Vector{Float64}
    down_offsets::Vector{Float64}
    heights::Vector{Float64}
end
_SpeciesSiteWorkspace(v) = _SpeciesSiteWorkspace(
    Int32[], Int32[], Int32[], Int32[], zeros(Int32, v), Int32[],
    _SpeciesChildren(Int32[], Int32[]), Float64[], Float64[], Float64[], Float64[])

function _species_site!(workspace, plan, tips)
    if plan.index === nothing || (plan.strategy !== :virtual &&
        !_species_site_virtual(plan.tree.nnodes, length(tips)))
        return plan.tree, plan.geometry, tips
    end
    index, tree, geometry = plan.index, plan.tree, plan.geometry
    selected, nodes = workspace.selected, workspace.nodes
    r = length(tips)
    resize!(selected, r)
    copyto!(selected, tips)
    sort!(selected; by = node -> index.position[node])
    resize!(nodes, 2r - 1)
    copyto!(nodes, 1, selected, 1, r)
    for i in 1:(r - 1)
        nodes[r + i] = _species_lca(index, tree, selected[i], selected[i + 1])
    end
    sort!(nodes; by = node -> index.position[node])
    unique!(nodes)
    m = length(nodes)
    for array in (workspace.parents, workspace.stack, workspace.children.first,
        workspace.children.next, workspace.factors, workspace.offsets, workspace.down_offsets)
        resize!(array, m)
    end
    resize!(workspace.heights, isempty(geometry.heights) ? 0 : m)
    fill!(workspace.children.first, 0)
    fill!(workspace.children.next, 0)
    top = 0
    for (i, node) in enumerate(nodes)
        workspace.local_index[node] = i
        while top > 0 && !_species_ancestor(index, nodes[workspace.stack[top]], node)
            top -= 1
        end
        parent = top == 0 ? Int32(0) : workspace.stack[top]
        workspace.parents[i] = parent
        f, up, down = parent == 0 ? _SPECIES_IDENTITY : _species_path(index, tree, nodes[parent], node)
        workspace.factors[i], workspace.offsets[i], workspace.down_offsets[i] = f, up, down
        !isempty(geometry.heights) && (workspace.heights[i] = geometry.heights[node])
        if parent != 0
            workspace.children.next[i] = workspace.children.first[parent]
            workspace.children.first[parent] = i
        end
        top += 1
        workspace.stack[top] = i
    end
    resize!(workspace.tips, r)
    for (i, node) in enumerate(tips)
        workspace.tips[i] = workspace.local_index[node]
    end
    topology = (; nnodes = m, root = Int32(1), parent_of_node = workspace.parents,
        children = workspace.children, postorder = m:-1:1, preorder = 1:m, tip_ids = Base.OneTo(m))
    local_geometry = (; factors = workspace.factors, offsets = workspace.offsets,
        down_offsets = workspace.down_offsets, heights = workspace.heights, scaled = geometry.scaled)
    return topology, local_geometry, workspace.tips
end

function _species_pse_values(rows, plan)
    if plan.index === nothing && plan.tree isa CompactTree
        return _pse_values(rows, plan.tree, plan.geometry, zeros(plan.tree.nnodes), zeros(plan.tree.nnodes))
    end
    result = fill(NaN, length(plan.rows))
    site_workspace = _SpeciesSiteWorkspace(plan.index === nothing ? 0 : plan.tree.nnodes)
    mass, means = Float64[], Float64[]
    precise_mass, precise_means = BigFloat[], BigFloat[]
    for (site, tips) in enumerate(plan.rows)
        length(tips) > 1 || continue
        tree, geometry, tips = _species_site!(site_workspace, plan, tips)
        if rows.precise[site] === nothing
            resize!(mass, tree.nnodes)
            resize!(means, tree.nnodes)
            result[site] = _pse_moments!(tree, tips, rows.weights[site], geometry, mass, means)
        else
            resize!(precise_mass, tree.nnodes)
            resize!(precise_means, tree.nnodes)
            result[site] = Float64(_pse_moments!(tree, tips, rows.precise[site], geometry, precise_mass, precise_means))
        end
    end
    return result
end

function _species_psc_values(plan)
    if plan.index === nothing && plan.tree isa CompactTree
        return _psc_values(plan.rows, plan.tree, plan.geometry)
    end
    result = fill(NaN, length(plan.rows))
    site_workspace = _SpeciesSiteWorkspace(plan.index === nothing ? 0 : plan.tree.nnodes)
    workspace = _ComdistntWorkspace(0)
    for (site, tips) in enumerate(plan.rows)
        length(tips) > 1 || continue
        tree, geometry, tips = _species_site!(site_workspace, plan, tips)
        for array in (workspace.first_distance, workspace.second_distance, workspace.first_tip, workspace.second_tip)
            resize!(array, tree.nnodes)
        end
        result[site] = _psc_value!(workspace, tree, tips, geometry)
    end
    return result
end

# A function boundary specializes all three kernels on either a full tree or
# the site's small topology. PSV and PSE reuse the same mean/mass storage;
# uniform PSE reuses PSV's value. PSC keeps its two distinct identities.
function _species_psd_site!(workspace, tree, tips, geometry, weights, precise, uniform)
    for array in (workspace.counts, workspace.means, workspace.mass,
        workspace.nearest.first_distance, workspace.nearest.second_distance,
        workspace.nearest.first_tip, workspace.nearest.second_tip)
        resize!(array, tree.nnodes)
    end
    value, _ = _psv_moments!(tree, tips, geometry, workspace.counts, workspace.means, workspace.mass)
    variability = geometry.scaled ? value : value / 2
    if uniform
        evenness = variability
    elseif precise === nothing
        evenness = _pse_moments!(tree, tips, weights, geometry, workspace.mass, workspace.means)
    else
        resize!(workspace.precise_mass, tree.nnodes)
        resize!(workspace.precise_means, tree.nnodes)
        evenness = Float64(_pse_moments!(tree, tips, precise, geometry,
            workspace.precise_mass, workspace.precise_means))
    end
    clustering = _psc_value!(workspace.nearest, tree, tips, geometry)
    return variability, evenness, clustering
end

function _species_psd_values(rows, plan, variance_by_richness)
    n = length(plan.rows)
    variability, evenness, clustering, variance = (fill(NaN, n) for _ in 1:4)
    richness = length.(rows.tips)
    site_workspace = _SpeciesSiteWorkspace(plan.index === nothing ? 0 : plan.tree.nnodes)
    workspace = (; counts = Int[], means = Float64[], mass = Float64[],
        precise_mass = BigFloat[], precise_means = BigFloat[], nearest = _ComdistntWorkspace(0))
    for (site, tips) in enumerate(plan.rows)
        length(tips) > 1 || continue
        tree, geometry, tips = _species_site!(site_workspace, plan, tips)
        v, e, c = _species_psd_site!(workspace, tree, tips, geometry,
            rows.weights[site], rows.precise[site], rows.uniform[site])
        variability[site], evenness[site], clustering[site] = v, e, c
        !isempty(variance_by_richness) && (variance[site] = variance_by_richness[richness[site]])
    end
    return (psv = variability, var_psv = variance, psc = clustering,
        psr = variability .* richness, var_psr = variance .* richness .^ 2,
        pse = evenness, richness)
end
