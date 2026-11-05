# Tree-dependent buffers survive only until the next tree. Community columns
# and pair schedules are fixed once; no source tree or per-tree beta is kept.
mutable struct _MeanBetaWorkspace
    branches::Vector{Vector{Int32}}
    sites::Vector{Int}
    edge_of_child::Vector{Int32}
    seen::Vector{Int}
    pointers::Vector{Int}
    cursors::Vector{Int}
    members::Vector{Int}
    shared::Matrix{Float64}
end

function _mean_beta_prepare(community, sample_labels, pairs)
    _validate_community_values(community)
    labels = _validated_labels(sample_labels, size(community, 1), "sample")
    selected = _community_pairs(pairs, length(labels))
    unique_pairs, outputs = Tuple{Int, Int}[], Int[]
    active = selected === nothing ? trues(length(labels)) : falses(length(labels))
    if selected !== nothing
        lookup = Dict{Tuple{Int, Int}, Int}()
        for (i, j) in selected
            key = minmax(i, j)
            index = get!(lookup, key) do
                push!(unique_pairs, key)
                length(unique_pairs)
            end
            push!(outputs, index)
            active[i] = active[j] = true
        end
    end
    rows = [Int32[] for _ in labels]
    _community_each_positive(community) do row, col, _
        active[row] && push!(rows[row], Int32(col))
    end
    workspace = _MeanBetaWorkspace([Int32[] for _ in labels], findall(active),
        Int32[], Int[], Int[], Int[], Int[], zeros(0, 0))
    return (; labels, selected, unique_pairs, outputs, rows, workspace)
end

function _mean_beta_branches!(workspace, tree, rows, species, complete)
    _validate_community_branch_lengths(tree)
    complete && length(species) != tree.ntips &&
        throw(ArgumentError("matrix columns must contain the complete tip set of every tree"))
    nodes = species == tree.tip_labels ? tree.tip_ids :
        Int32[get(tree.tipname_to_id, label, 0) for label in species]
    all(node -> node != 0 && tree.is_tip[node], nodes) ||
        throw(ArgumentError("community species must be present in every tree"))
    resize!(workspace.edge_of_child, tree.nnodes)
    resize!(workspace.seen, tree.nnodes)
    fill!(workspace.seen, 0)
    for edge in 1:tree.nedges
        workspace.edge_of_child[tree.child_of_edge[edge]] = edge
    end
    for site in workspace.sites
        edges = workspace.branches[site]
        empty!(edges)
        for col in rows[site]
            node = nodes[col]
            while node != tree.root && workspace.seen[node] != site
                workspace.seen[node] = site
                push!(edges, workspace.edge_of_child[node])
                node = tree.parent_of_node[node]
            end
        end
        sort!(edges)
    end
    return workspace.branches
end

# Invert the sparse site/branch membership directly. Site indices in each
# edge's range are sorted, so the lower-triangle scatter can stop early.
# Only indices are needed; retain their capacities across posterior trees
# instead of building two generic CSC matrices with unused numeric values.
function _mean_beta_members!(workspace, nedges)
    pointers, cursors, members = workspace.pointers, workspace.cursors, workspace.members
    resize!(pointers, nedges + 1)
    fill!(pointers, 0)
    pointers[1] = 1
    for edges in workspace.branches, edge in edges
        pointers[edge + 1] += 1
    end
    for edge in 1:nedges
        pointers[edge + 1] += pointers[edge]
    end
    resize!(cursors, nedges)
    copyto!(cursors, 1, pointers, 1, nedges)
    resize!(members, pointers[end] - 1)
    for (site, edges) in enumerate(workspace.branches), edge in edges
        members[cursors[edge]] = site
        cursors[edge] += 1
    end
    return pointers, members
end

# The same scaled formulas as _community_beta_components, evaluating only
# the selected component. Nestedness needs both ratios and the same NaN mask.
@inline function _mean_beta_component(a, b, c, family, component)
    component == 3 && return _community_beta_components(a, b, c, family)[3]
    if component == 1
        a, b, c = _community_scaled_components(a, b, c)
        denominator = family === :sorensen ? 2a + b + c : a + b + c
        return denominator > 0 ? (b + c) / denominator : NaN
    end
    unique = min(b, c)
    scale = max(a, unique)
    a, unique = scale > 0 ? (a / scale, unique / scale) : (0.0, 0.0)
    denominator = family === :sorensen ? a + unique : a + 2unique
    return denominator > 0 ? (family === :sorensen ? unique : 2unique) / denominator : NaN
end

function _mean_beta_visit!(state, workspace, lengths, selected, family, component, policy, block_size)
    branches = workspace.branches
    n = length(branches)
    if selected !== nothing
        for (index, (lo, hi)) in enumerate(selected)
            abc = _community_branch_components(branches[hi], branches[lo], lengths)
            _mean_dist_add!(state, index, _mean_beta_component(abc..., family, component), policy)
        end
    elseif n < 128
        for j in 1:n, i in j:n
            abc = _community_branch_components(branches[i], branches[j], lengths)
            _mean_dist_add!(state, i + (j - 1) * n, _mean_beta_component(abc..., family, component), policy)
        end
    else
        # This is the same sparse incidence product, scattered directly into
        # one dense column block. No generic CSC/weighted copy, sparse result,
        # or per-block slicing/conversion is needed. Only the lower
        # triangle contributes; each cell adds original edges in sorted order.
        pointers, members = _mean_beta_members!(workspace, length(lengths))
        pd_values = zeros(n)
        for site in 1:n, edge in branches[site]
            pd_values[site] += lengths[edge]
        end
        size(workspace.shared) == (n, min(n, block_size)) ||
            (workspace.shared = zeros(n, min(n, block_size)))
        shared = workspace.shared
        for first in 1:block_size:n
            last = min(n, first + block_size - 1)
            fill!(shared, 0)
            @inbounds for j in first:last, edge in branches[j]
                edge_length = lengths[edge]
                iszero(edge_length) && continue
                for entry in (pointers[edge + 1] - 1):-1:pointers[edge]
                    i = members[entry]
                    i < j && break
                    shared[i, j - first + 1] += edge_length
                end
            end
            for j in first:last, i in j:n
                abc = _community_components_from_shared(branches, lengths,
                    shared[i, j - first + 1], pd_values, i, j)
                _mean_dist_add!(state, i + (j - 1) * n,
                    _mean_beta_component(abc..., family, component), policy)
            end
        end
    end
    return nothing
end

function _mean_beta(community, trees, species, sample_labels, complete;
    index_family::Symbol = :sorensen, component::Symbol = :total,
    pairs = nothing, nan_policy::Symbol = :skip, block_size = 256,
)
    _validate_community_beta_family(index_family)
    component_index = findfirst(==(component), (:total, :turnover, :nestedness))
    component_index === nothing && throw(ArgumentError("component must be :total, :turnover or :nestedness"))
    nan_policy in (:skip, :propagate) || throw(ArgumentError("nan_policy must be :skip or :propagate"))
    block_size isa Integer && !(block_size isa Bool) && 1 <= block_size <= typemax(Int) ||
        throw(ArgumentError("block_size must be a positive integer"))
    plan = _mean_beta_prepare(community, sample_labels, pairs)
    n = length(plan.labels)
    selected = plan.selected === nothing ? nothing : plan.unique_pairs
    shape = selected === nothing ? (n, n) : (length(selected),)
    state = nothing
    replicates = 0
    for tree in trees
        tree isa CompactTree || throw(ArgumentError("trees must yield CompactTree values"))
        if species === nothing
            size(community, 2) == tree.ntips ||
                throw(ArgumentError("matrix must have one column per tip in the first tree"))
            species = _validated_labels(tree.tip_labels, tree.ntips, "tree tip")
        end
        branches = _mean_beta_branches!(plan.workspace, tree, plan.rows, species, complete)
        if state === nothing
            state = (; sums = zeros(shape), counts = zeros(Int, shape), scales = zeros(shape), corrections = zeros(shape))
        end
        _mean_beta_visit!(state, plan.workspace, tree.edge_length, selected, index_family,
            component_index, nan_policy, min(Int(block_size), max(n, 1)))
        replicates += 1
    end
    replicates > 0 || throw(ArgumentError("trees must not be empty"))
    _mean_dist_finish!(state, selected, n)
    distances = selected === nothing ? state.sums : state.sums[plan.outputs]
    counts = selected === nothing ? state.counts : state.counts[plan.outputs]
    return MeanDistResult(plan.labels, distances, counts, plan.selected, replicates, component, index_family)
end

"""
    mean_dist(community, trees; component=:total, index_family=:sorensen,
              pairs=nothing, nan_policy=:skip, block_size=256,
              tip_labels=nothing, sample_labels=nothing)
    mean_dist(data::CommunityData, trees; kwargs...)

Compute one phylogenetic beta component over a nonempty, single-pass tree
iterator and accumulate directly into `MeanDistResult`. No per-tree beta
result matrices/vectors are retained. Full output uses blocked shared-branch
sparse multiplication; `block_size` limits its temporary columns. Selected
pairs evaluate each unordered request once and retain only requested endpoints.

Matrix columns must describe the complete tip set of every tree. `tip_labels`
fixes their meaning; when omitted, the first tree's labels define it for all
later trees, even when tip order changes. `CommunityData` supplies both label
sets and may omit unobserved tree species. Unknown species fail, including
zero columns. Every tree is validated and matched separately; no pruning or
reuse of another tree's topology occurs. Presence records are prepared once.

Labels/pairs, component-specific NaNs, `nan_policy`, finite counts and result
ownership follow `mean_dist(inputs)`. Self pairs are zero or NaN. Use an
explicit `read_trees(...) do trees ... end` scope for file streams. Complete
output and four accumulator buffers still need O(sites²) space; selected
output uses O(sites + requests) in addition to the current tree/branch sets.
"""
function mean_dist(community::AbstractMatrix{<:Real}, trees;
    tip_labels = nothing, sample_labels = nothing, kwargs...,
)
    species = tip_labels === nothing ? nothing : _validated_labels(tip_labels, size(community, 2), "tip")
    return _mean_beta(community, trees, species, sample_labels, true; kwargs...)
end

function mean_dist(data::CommunityData, trees; kwargs...)
    _community_data_keywords(kwargs)
    species = _validated_labels(data.species_labels, size(data, 2), "species")
    return _mean_beta(data.matrix, trees, species, data.sample_labels, false; kwargs...)
end
