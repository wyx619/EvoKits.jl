# SES shares one sparse null-model engine; only a single replicate is retained.
# These workspaces scale with observed occurrences and tree size, not the
# product of site and species counts.

# Shared finite-draw Welford statistics for site and pairwise SES. Rank counts
# start at one for the observed value; callers choose the documented p-value
# denominator (historical runs+1 for site SES, finite count+1 for beta SES).
struct _SESOnlineStats{N}
    counts::Array{Int, N}
    means::Array{Float64, N}
    m2::Array{Float64, N}
    ranks::Array{Float64, N}
end

_SESOnlineStats(shape::Tuple) = _SESOnlineStats(zeros(Int, shape),
    zeros(Float64, shape), zeros(Float64, shape), ones(Float64, shape))

@inline function _ses_online_update!(stats, i, observed, value)
    isfinite(value) || return nothing
    stats.counts[i] += 1
    delta = value - stats.means[i]
    stats.means[i] += delta / stats.counts[i]
    stats.m2[i] += delta * (value - stats.means[i])
    if isfinite(observed)
        stats.ranks[i] += value < observed ? 1.0 : value == observed ? 0.5 : 0.0
    end
    return nothing
end

function _ses_online_finish!(stats, observed, indices = eachindex(observed))
    sds, z = (fill(NaN, size(observed)) for _ in 1:2)
    for i in indices
        stats.counts[i] == 0 && (stats.means[i] = NaN)
        sds[i] = stats.counts[i] > 1 ? sqrt(max(0.0, stats.m2[i] / (stats.counts[i] - 1))) : NaN
        z[i] = isfinite(sds[i]) && sds[i] > 0 ? (observed[i] - stats.means[i]) / sds[i] : NaN
    end
    return (; means = stats.means, sds, z, ranks = stats.ranks, counts = stats.counts)
end

mutable struct _SesPDWorkspace
    edge_of_child::Vector{Int32}
    seen::Vector{Int}
    generation::Int
    index::Union{Nothing, _TreeQueryIndex}
    include_root::Bool
end

function _SesPDWorkspace(tree, include_root)
    return _SesPDWorkspace(_community_edge_index(tree), zeros(Int, tree.nnodes), 0,
        include_root ? nothing : _build_tree_query_index(tree), include_root)
end

function _ses_pd_site!(tree, tips, workspace)
    isempty(tips) && return 0.0
    !workspace.include_root && length(tips) == 1 && return NaN
    if workspace.generation == typemax(Int)
        fill!(workspace.seen, 0)
        workspace.generation = 0
    end
    workspace.generation += 1
    root = tree.root
    if !workspace.include_root
        root = tree.tip_ids[first(tips)]
        for tip in Iterators.drop(tips, 1)
            root = _tree_lca(workspace.index, root, tree.tip_ids[tip])
        end
    end
    value = 0.0
    for tip in tips
        node = tree.tip_ids[tip]
        while node != root && workspace.seen[node] != workspace.generation
            workspace.seen[node] = workspace.generation
            value += tree.edge_length[workspace.edge_of_child[node]]
            node = tree.parent_of_node[node]
        end
    end
    return value
end

struct _SesDistanceWorkspace
    index::_TreeQueryIndex
    virtual::_NriMetricsWorkspace
    edge_of_child::Vector{Int32}
    node_weights::Vector{Float64}
    masses::Vector{Float64}
    weights::Vector{Float64}
end

function _SesDistanceWorkspace(tree, weighted)
    return _SesDistanceWorkspace(_build_tree_query_index(tree), _NriMetricsWorkspace(tree), _community_edge_index(tree),
        weighted ? zeros(tree.nnodes) : Float64[], Float64[], Float64[])
end

function _ses_abundance_rows(community, perm, rows)
    inverse = invperm(perm)
    weights = [Float64[community[site, inverse[tip]] for tip in row] for (site, row) in enumerate(rows)]
    all(row -> all(isfinite, row), weights) ||
        throw(ArgumentError("SES abundance weights must be representable as finite Float64 values"))
    return weights
end

# Preserve the existing abundance null-model convention: randomize presence,
# keep original weights at their tree-tip positions. Null tips absent from the
# original site have zero weight, but remain nearest-neighbor candidates.
function _ses_site_weights!(dest, tips, original, weights)
    resize!(dest, length(tips))
    j = 1
    for (i, tip) in enumerate(tips)
        while j <= length(original) && original[j] < tip
            j += 1
        end
        dest[i] = j <= length(original) && original[j] == tip ? weights[j] : 0.0
    end
    return dest
end

function _ses_distance_site!(tree, tips, workspace, metric, weights)
    isempty(tips) && return NaN
    if length(tips) == 1
        return metric === :mpd && weights !== nothing && first(weights) > 0 ? 0.0 : NaN
    end
    mpd, mntd = _nri_community_metrics!(tree, tips, workspace.index, workspace.virtual, workspace.edge_of_child)
    weights === nothing && return metric === :mpd ? mpd : mntd
    total = sum(weights)
    total > 0 || return NaN
    for (tip, weight) in zip(tips, weights)
        workspace.node_weights[tree.tip_ids[tip]] = weight
    end
    virtual = workspace.virtual
    nodes = virtual.virtual_nodes
    mass = workspace.masses
    resize!(mass, length(nodes))
    fill!(mass, 0.0)
    nearest = 0.0
    for (i, node) in enumerate(nodes)
        if virtual.marks[node] == virtual.generation
            weight = workspace.node_weights[node]
            mass[i] = weight
            nearest += weight * virtual.updist[i]
        end
    end
    metric === :mntd && return nearest / total
    pair_sum = 0.0
    for i in length(nodes):-1:2
        parent = Int(virtual.parent_index[i])
        edge_length = virtual.edge_lengths[i]
        pair_sum += edge_length * mass[i] * (total - mass[i])
        mass[parent] += mass[i]
    end
    return 2pair_sum / (total * total)
end

function _ses_evaluate!(values, metric, tree, rows, workspace::_SesPDWorkspace, observed, weights)
    for (site, row) in enumerate(rows)
        values[site] = _ses_pd_site!(tree, row, workspace)
    end
    return values
end

function _ses_evaluate!(values, metric, tree, rows, workspace::_SesDistanceWorkspace, observed, weights)
    for (site, row) in enumerate(rows)
        site_weights = weights === nothing ? nothing :
            _ses_site_weights!(workspace.weights, row, observed[site], weights[site])
        values[site] = _ses_distance_site!(tree, row, workspace, metric, site_weights)
    end
    return values
end

function _ses_summary(metric, community, tree;
    null_model, runs, iterations, include_root = true, abundance_weighted = false,
    tip_labels, sample_labels, rng, max_attempts_per_swap,
)
    runs >= 1 || throw(ArgumentError("runs must be positive"))
    iterations >= 0 || throw(ArgumentError("iterations must be non-negative"))
    max_attempts_per_swap >= 1 || throw(ArgumentError("max_attempts_per_swap must be positive"))
    model = _tools_normalize_phylo_null(null_model)
    _validate_community_branch_lengths(tree)
    size(community, 2) == tree.ntips || throw(ArgumentError("community matrix must have one column per tree tip"))
    labels = _tools_sample_labels(size(community, 1), sample_labels)
    perm = _tools_tip_permutation(tree, tip_labels)
    observed = _community_tip_rows(community, perm)
    weights = abundance_weighted ? _ses_abundance_rows(community, perm, observed) : nothing
    null = _CommunityNullWorkspace(observed, tree.ntips)
    workspace = metric === :pd ? _SesPDWorkspace(tree, include_root) :
        _SesDistanceWorkspace(tree, abundance_weighted)
    obs, random = (Vector{Float64}(undef, length(observed)) for _ in 1:2)
    _ses_evaluate!(obs, metric, tree, observed, workspace, observed, weights)
    stats = _SESOnlineStats(size(obs))
    for _ in 1:runs
        rows = _randomize_community!(null, observed, model, iterations, max_attempts_per_swap, rng)
        _ses_evaluate!(random, metric, tree, rows, workspace, observed, weights)
        for i in eachindex(random)
            _ses_online_update!(stats, i, obs[i], random[i])
        end
    end
    (; means, sds, ranks, z) = _ses_online_finish!(stats, obs)
    return (; labels, richness = null.richness, obs, means, sds, ranks, z,
        p = ranks ./ (runs + 1), model)
end
