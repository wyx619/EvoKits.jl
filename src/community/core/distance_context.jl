"""
    CommunityDistanceContext

Owned community/tree snapshot for repeated `comdist` and `comdistnt` calls,
created with [`prepare_community_distances`](@ref). Stores presence and abundance
records, one full or compressed topology, and reusable workspaces. Its fields
are internal and must not be mutated. Calls sharing a context are serialized;
returned arrays are independent. No distance matrix is cached.
"""
struct CommunityDistanceContext{R, U, T, P}
    sample_labels::Vector{String}
    rows::R
    presence::U
    topology::T
    pairs::Union{Nothing, Vector{Tuple{Int, Int}}}
    order::Vector{Int}
    plan::P
    mean_workspace::_ComdistWorkspace
    nearest_workspace::_ComdistntWorkspace
    lock::ReentrantLock
end

# Uniform and abundance paths share node IDs and branch lengths, but keep
# their original raw weights and direction-total definitions independently.
function _distance_presence_rows(rows)
    raw = [ones(length(nodes)) for nodes in rows.nodes]
    weights = [isempty(row) ? Float64[] : fill(1 / length(row), length(row)) for row in raw]
    scales = [isempty(row) ? 0.0 : 1.0 for row in raw]
    totals = Float64.(length.(raw))
    return merge(rows, (; raw, weights, scales, totals))
end

"""
    prepare_community_distances(community, tree; tip_labels=nothing,
                                sample_labels=nothing, pairs=nothing)
    prepare_community_distances(data::CommunityData, tree; pairs=nothing)

Validate and snapshot a community collection and tree once for repeated
`comdist(context; abundance_weighted=false, pairs=...)` and
`comdistnt(context; abundance_weighted=false, exclude_conspecifics=false, pairs=...)`.
Both weighting modes share prepared labels, node records, branch lengths and
topology. Inputs may be dense/CSC matrices or `CommunityData`; matrix columns
must cover all tree tips. CommunityData aligns its species labels as usual.
Site labels must be unique and nonempty; defaults are `sample_1`, etc.

Without preparation-time `pairs`, later calls may request complete matrices
or any subset of pairs using the topology for all prepared sites. Supplying
`pairs` here prepares just their nonempty endpoints, fixes the request order
and caches both schedules. Later calls then default to that vector request;
changing it requires rebuilding the context. Original site indices apply,
including repeated, reversed and self pairs. Empty requests are supported.

All input values, labels and branch lengths are checked even for unrequested
sites. Later changes to input values, labels or tree topology do not affect
the snapshot; rebuild it for changed data. The original matrix/tree is not
retained. Records and workspaces need linear auxiliary storage in sites,
retained observations and working-tree nodes. Full outputs still need O(n²)
space; selected outputs use O(k) space. Each returned array is newly owned.

This context is separate from the presence-only `CommunityContext`, which
does not retain abundances or distance topology. Reusing a context skips
input validation, alignment and tree construction; only new pair arguments
are validated. Fields are internal and must not be mutated. Calls sharing
its mutable workspaces are protected by a lock and run one at a time.
"""
function prepare_community_distances(community::AbstractMatrix{<:Real}, tree::CompactTree;
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    pairs = nothing,
)
    labels = _validated_labels(sample_labels, size(community, 1), "sample")
    selected = _community_pairs(pairs, size(community, 1))
    active = selected === nothing ? nothing : falses(size(community, 1))
    if selected !== nothing
        for (i, j) in selected
            active[i] = active[j] = true
        end
    end
    rows = _community_distance_rows(community, tree, tip_labels, true;
        active_sites = active, need_lengths = false)
    # BigInt/BigFloat and abstract Real containers may hold mutable numbers.
    # Own those values as well as their vectors, without copying bitstype rows.
    if !isbitstype(eltype(rows.scales))
        raw, scales = deepcopy(rows.raw), copy(rows.scales)
        for i in eachindex(raw)
            scales[i] = isempty(raw[i]) ? zero(eltype(scales)) : maximum(raw[i])
        end
        rows = merge(rows, (; raw, scales))
    end
    order = selected === nothing ? Int[] : sortperm(selected; by = pair -> minmax(pair...))
    plan = selected === nothing ? nothing : _comdistnt_pair_plan(rows, selected)
    if selected !== nothing
        for i in eachindex(rows.nodes)
            plan.used[i] && continue
            # Replace populated vectors rather than empty! them: a dense
            # endpoint with an empty partner must not retain its old capacity.
            if !isempty(rows.nodes[i])
                rows.nodes[i] = Int32[]
                rows.raw[i] = empty(rows.raw[i])
                rows.weights[i] = Float64[]
            end
            rows.scales[i] = zero(eltype(rows.scales))
            rows.totals[i] = 0
        end
    end
    targets = count(!isempty, rows.nodes)
    if targets == 0 || _comdist_use_virtual(rows, tree, targets)
        rows, topology = _comdist_virtual(rows, tree)
    else
        rows = merge(rows, (; lengths = _comdist_lengths(tree)))
        topology = (; nnodes = tree.nnodes, root = tree.root,
            parent_of_node = copy(tree.parent_of_node), children = deepcopy(tree.children),
            postorder = copy(tree.postorder), preorder = copy(tree.preorder))
    end
    presence = _distance_presence_rows(rows)
    return CommunityDistanceContext(labels, rows, presence, topology, selected, order, plan,
        _ComdistWorkspace(topology.nnodes), _ComdistntWorkspace(topology.nnodes), ReentrantLock())
end

function _distance_context_pairs(context, pairs)
    pairs === context.pairs && return pairs
    selected = _community_pairs(pairs, length(context.sample_labels))
    if context.pairs !== nothing
        selected == context.pairs ||
            throw(ArgumentError("context has fixed pairs; rebuild without preparation-time pairs to change requests"))
        return context.pairs
    end
    return selected
end

"""`comdist(context::CommunityDistanceContext; abundance_weighted=false, pairs=...)` reuses prepared records, topology and workspace. Fixed contexts default to their prepared pairs."""
function comdist(context::CommunityDistanceContext; abundance_weighted::Bool = false,
    pairs = context.pairs)
    selected = _distance_context_pairs(context, pairs)
    rows = abundance_weighted ? context.rows : context.presence
    order = selected === context.pairs ? context.order : sortperm(selected; by = pair -> minmax(pair...))
    lock(context.lock)
    try
        return selected === nothing ? _comdist(rows, context.topology, rows.weights, context.mean_workspace) :
            _comdist_selected(rows, context.topology, selected, order, context.mean_workspace)
    finally
        unlock(context.lock)
    end
end

"""`comdistnt(context::CommunityDistanceContext; abundance_weighted=false, exclude_conspecifics=false, pairs=...)` reuses prepared records, topology and workspace. Fixed contexts default to their prepared pairs."""
function comdistnt(context::CommunityDistanceContext; abundance_weighted::Bool = false,
    exclude_conspecifics::Bool = false, pairs = context.pairs)
    selected = _distance_context_pairs(context, pairs)
    rows = abundance_weighted ? context.rows : context.presence
    plan = selected === context.pairs ? context.plan : _comdistnt_pair_plan(rows, selected)
    lock(context.lock)
    try
        return selected === nothing ?
            _comdistnt(rows, context.topology, abundance_weighted, exclude_conspecifics, context.nearest_workspace) :
            _comdistnt_selected(rows, context.topology, plan, abundance_weighted, exclude_conspecifics, context.nearest_workspace)
    finally
        unlock(context.lock)
    end
end
