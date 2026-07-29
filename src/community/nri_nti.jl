"""
    NriNtiResult

NRI/NTI result for one observed community and its local null-model pool.
"""
Base.@kwdef struct NriNtiResult
    observed_mpd::Float64 = NaN
    observed_mntd::Float64 = NaN
    null_mpd_mean::Float64 = NaN
    null_mpd_sd::Float64 = NaN
    null_mntd_mean::Float64 = NaN
    null_mntd_sd::Float64 = NaN
    nri::Float64 = NaN
    nti::Float64 = NaN
    pool_richness::Int = 0
    observed_richness::Int = 0
    nreps::Int = 0
    seed::UInt64 = UInt64(0)
    status::Symbol = :ok
end

function Base.getproperty(result::NriNtiResult, name::Symbol)
    name === :MPD_obs && return getfield(result, :observed_mpd)
    name === :MNTD_obs && return getfield(result, :observed_mntd)
    name === :MPD_null_mean && return getfield(result, :null_mpd_mean)
    name === :MPD_null_sd && return getfield(result, :null_mpd_sd)
    name === :MNTD_null_mean && return getfield(result, :null_mntd_mean)
    name === :MNTD_null_sd && return getfield(result, :null_mntd_sd)
    name === :NRI && return getfield(result, :nri)
    name === :NTI && return getfield(result, :nti)
    return getfield(result, name)
end

function _nri_tip_positions(
    tree::CompactTree,
    values,
    name::AbstractString;
    tip_position_of_node::Union{Nothing,AbstractVector{<:Integer}} = nothing,
)
    output = if eltype(values) <: Integer
        Int.(values)
    else
        positions = Int[]
        for label0 in values
            label = String(label0)
            node = get(tree.tipname_to_id, label, Int32(0))
            node == 0 && throw(ArgumentError("$name contains unknown tree tip: $label"))
            position = tip_position_of_node === nothing ?
                findfirst(==(node), tree.tip_ids) : tip_position_of_node[Int(node)]
            (position === nothing || position == 0) && throw(ArgumentError("$name tip is not indexed: $label"))
            push!(positions, Int(position))
        end
        positions
    end
    all((1 .<= output) .& (output .<= tree.ntips)) || throw(ArgumentError("$name contains tip positions outside 1:$(tree.ntips)"))
    length(unique(output)) == length(output) || throw(ArgumentError("$name contains duplicate tips"))
    return output
end

"""Reusable virtual-tree buffers for the repeated NRI/NTI null samples."""
mutable struct _NriMetricsWorkspace
    selected_nodes::Vector{Int32}
    virtual_nodes::Vector{Int32}
    parent_index::Vector{Int32}
    stack::Vector{Int32}
    first_child::Vector{Int32}
    next_sibling::Vector{Int32}
    counts::Vector{Int32}
    down::Vector{Float64}
    updist::Vector{Float64}
    marks::Vector{Int32}
    generation::Int32
end

function _NriMetricsWorkspace(tree::CompactTree)
    return _NriMetricsWorkspace(
        Int32[], Int32[], Int32[], Int32[], Int32[], Int32[], Int32[],
        Float64[], Float64[], zeros(Int32, tree.nnodes), Int32(0),
    )
end

function _nri_next_generation!(workspace::_NriMetricsWorkspace)
    if workspace.generation == typemax(Int32)
        fill!(workspace.marks, Int32(0))
        workspace.generation = Int32(1)
    else
        workspace.generation += Int32(1)
    end
    return workspace.generation
end

"""Compute non-abundance MPD and MNTD without per-call container allocation."""
function _nri_community_metrics!(
    tree::CompactTree,
    tip_positions::AbstractVector{<:Integer},
    index::_TreeQueryIndex,
    workspace::_NriMetricsWorkspace,
)
    richness = length(tip_positions)
    selected_nodes = workspace.selected_nodes
    resize!(selected_nodes, richness)
    @inbounds for i in 1:richness
        selected_nodes[i] = tree.tip_ids[Int(tip_positions[i])]
    end
    sort!(selected_nodes; by = node -> index.tin[Int(node)])

    virtual_nodes = workspace.virtual_nodes
    resize!(virtual_nodes, 2 * richness - 1)
    copyto!(virtual_nodes, 1, selected_nodes, 1, richness)
    @inbounds for i in 1:(richness - 1)
        virtual_nodes[richness + i] = _tree_lca(index, selected_nodes[i], selected_nodes[i + 1])
    end
    sort!(virtual_nodes; by = node -> index.tin[Int(node)])
    unique!(virtual_nodes)
    m = length(virtual_nodes)

    generation = _nri_next_generation!(workspace)
    @inbounds for node in selected_nodes
        workspace.marks[Int(node)] = generation
    end

    parent_index = workspace.parent_index
    stack = workspace.stack
    resize!(parent_index, m)
    resize!(stack, m)
    stack_length = 0
    @inbounds for i in 1:m
        node = virtual_nodes[i]
        while stack_length > 0 && !_tree_is_ancestor(index, virtual_nodes[Int(stack[stack_length])], node)
            stack_length -= 1
        end
        parent_index[i] = stack_length == 0 ? Int32(0) : stack[stack_length]
        stack_length += 1
        stack[stack_length] = Int32(i)
    end

    first_child = workspace.first_child
    next_sibling = workspace.next_sibling
    resize!(first_child, m)
    resize!(next_sibling, m)
    fill!(first_child, Int32(0))
    fill!(next_sibling, Int32(0))
    @inbounds for i in 2:m
        parent = Int(parent_index[i])
        next_sibling[i] = first_child[parent]
        first_child[parent] = Int32(i)
    end

    counts = workspace.counts
    down = workspace.down
    updist = workspace.updist
    resize!(counts, m)
    resize!(down, m)
    resize!(updist, m)
    fill!(counts, Int32(0))
    fill!(down, Inf)
    fill!(updist, Inf)
    @inbounds for i in 1:m
        if workspace.marks[Int(virtual_nodes[i])] == generation
            counts[i] = Int32(1)
            down[i] = 0.0
        end
    end

    pair_sum = 0.0
    @inbounds for i in m:-1:2
        parent = Int(parent_index[i])
        edge_length = tree.dist_from_root[Int(virtual_nodes[i])] - tree.dist_from_root[Int(virtual_nodes[parent])]
        count = counts[i]
        pair_sum += edge_length * count * (richness - count)
        counts[parent] += count
        down[parent] = min(down[parent], down[i] + edge_length)
    end

    @inbounds for parent in 1:m
        best1 = workspace.marks[Int(virtual_nodes[parent])] == generation ? 0.0 : Inf
        best2 = Inf
        best_child = Int32(0)
        child = first_child[parent]
        while child != 0
            c = Int(child)
            edge_length = tree.dist_from_root[Int(virtual_nodes[c])] - tree.dist_from_root[Int(virtual_nodes[parent])]
            candidate = down[c] + edge_length
            if candidate < best1
                best2, best1, best_child = best1, candidate, child
            elseif candidate < best2
                best2 = candidate
            end
            child = next_sibling[c]
        end
        child = first_child[parent]
        while child != 0
            c = Int(child)
            edge_length = tree.dist_from_root[Int(virtual_nodes[c])] - tree.dist_from_root[Int(virtual_nodes[parent])]
            sibling_best = child == best_child ? best2 : best1
            updist[c] = edge_length + min(updist[parent], sibling_best)
            child = next_sibling[c]
        end
    end

    nearest_sum = 0.0
    @inbounds for i in 1:m
        workspace.marks[Int(virtual_nodes[i])] == generation || continue
        child_best = Inf
        child = first_child[i]
        while child != 0
            c = Int(child)
            edge_length = tree.dist_from_root[Int(virtual_nodes[c])] - tree.dist_from_root[Int(virtual_nodes[i])]
            child_best = min(child_best, down[c] + edge_length)
            child = next_sibling[c]
        end
        nearest_sum += min(updist[i], child_best)
    end

    return pair_sum / (richness * (richness - 1) / 2), nearest_sum / richness
end

mutable struct _NriOnlineMoments
    n::Int
    mean::Float64
    m2::Float64
end

_NriOnlineMoments() = _NriOnlineMoments(0, 0.0, 0.0)

function _nri_update!(moments::_NriOnlineMoments, value::Float64)
    moments.n += 1
    delta = value - moments.mean
    moments.mean += delta / moments.n
    moments.m2 += delta * (value - moments.mean)
    return moments
end

function _nri_sd(moments::_NriOnlineMoments)
    moments.n >= 2 || return NaN
    return sqrt(max(0.0, moments.m2 / (moments.n - 1)))
end

function _nri_nti_positions(
    tree::CompactTree,
    pool_positions::AbstractVector{<:Integer},
    observed_positions::AbstractVector{<:Integer},
    query_index::_TreeQueryIndex,
    workspace::_NriMetricsWorkspace;
    nreps::Integer,
    random::AbstractRNG,
    effective_seed::UInt64,
)
    observed_mpd, observed_mntd = _nri_community_metrics!(tree, observed_positions, query_index, workspace)
    mpd_moments = _NriOnlineMoments()
    mntd_moments = _NriOnlineMoments()
    working_pool = Int.(pool_positions)
    sampled = Vector{Int}(undef, length(observed_positions))
    swap_positions = Vector{Int}(undef, length(observed_positions))

    for _ in 1:nreps
        for i in eachindex(sampled)
            j = rand(random, i:length(working_pool))
            swap_positions[i] = j
            working_pool[i], working_pool[j] = working_pool[j], working_pool[i]
            sampled[i] = working_pool[i]
        end
        null_mpd, null_mntd = _nri_community_metrics!(tree, sampled, query_index, workspace)
        _nri_update!(mpd_moments, null_mpd)
        _nri_update!(mntd_moments, null_mntd)
        for i in length(sampled):-1:1
            j = swap_positions[i]
            working_pool[i], working_pool[j] = working_pool[j], working_pool[i]
        end
    end

    mpd_sd = _nri_sd(mpd_moments)
    mntd_sd = _nri_sd(mntd_moments)
    status = mpd_sd > 0.0 && mntd_sd > 0.0 ? :ok : :null_sd_zero
    nri = mpd_sd > 0.0 ? (mpd_moments.mean - observed_mpd) / mpd_sd : NaN
    nti = mntd_sd > 0.0 ? (mntd_moments.mean - observed_mntd) / mntd_sd : NaN
    return NriNtiResult(
        observed_mpd = observed_mpd,
        observed_mntd = observed_mntd,
        null_mpd_mean = mpd_moments.mean,
        null_mpd_sd = mpd_sd,
        null_mntd_mean = mntd_moments.mean,
        null_mntd_sd = mntd_sd,
        nri = nri,
        nti = nti,
        pool_richness = length(pool_positions),
        observed_richness = length(observed_positions),
        nreps = Int(nreps),
        seed = effective_seed,
        status = status,
    )
end

"""
    nri_nti(tree, pool, observed; nreps = 999, seed = 20260713, grid_id = nothing)

Compute NRI and NTI for an observed community sampled without replacement from
its local `pool`. `pool` and `observed` may contain tip positions or tip labels.
The observed community must contain at least two tips and be a subset of the
pool. Results are deterministic for the same `seed`. Pass `grid_id` for a
batch analysis to derive the same stable per-grid seed used by the Swq
validation workflow. `grid_id` cannot be combined with an injected `rng`.
"""
function nri_nti(
    tree::CompactTree,
    pool,
    observed;
    nreps::Integer = 999,
    seed::Integer = 20260713,
    rng::Union{Nothing, AbstractRNG} = nothing,
    grid_id::Union{Nothing, AbstractString} = nothing,
)
    nreps >= 2 || throw(ArgumentError("nreps must be at least 2 to calculate sample standard deviations"))
    rng !== nothing && grid_id !== nothing &&
        throw(ArgumentError("grid_id cannot be combined with an injected rng"))
    all(isfinite, tree.edge_length) && all(tree.edge_length .>= 0.0) || throw(ArgumentError("tree contains invalid branch lengths"))
    pool_values = collect(pool)
    observed_values = collect(observed)
    tip_position_of_node = if eltype(pool_values) <: Integer && eltype(observed_values) <: Integer
        nothing
    else
        index = zeros(Int, tree.nnodes)
        for (position, node) in enumerate(tree.tip_ids)
            index[Int(node)] = position
        end
        index
    end
    pool_positions = _nri_tip_positions(tree, pool_values, "pool"; tip_position_of_node)
    observed_positions = _nri_tip_positions(tree, observed_values, "observed"; tip_position_of_node)
    length(observed_positions) >= 2 || throw(ArgumentError("observed must contain at least two tips"))
    length(pool_positions) >= length(observed_positions) || throw(ArgumentError("pool is smaller than observed"))
    pool_set = Set(pool_positions)
    all(in(pool_set), observed_positions) || throw(ArgumentError("observed is not a subset of pool"))

    query_index = _build_tree_query_index(tree)
    workspace = _NriMetricsWorkspace(tree)
    effective_seed = grid_id === nothing ? UInt64(seed) : stable_grid_seed(seed, grid_id)
    random = rng === nothing ? StableRNG(grid_id === nothing ? seed : effective_seed) : rng
    return _nri_nti_positions(tree, pool_positions, observed_positions, query_index, workspace;
                              nreps, random, effective_seed)
end
