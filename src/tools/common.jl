@inline function _tools_normalize_phylo_null(null_model::Symbol)
    null_model === :taxa_labels && return null_model
    null_model === Symbol("taxa.labels") && return :taxa_labels
    null_model === :richness && return null_model
    null_model === :frequency && return null_model
    null_model === :sample_pool && return null_model
    null_model === Symbol("sample.pool") && return :sample_pool
    null_model === :phylogeny_pool && return null_model
    null_model === Symbol("phylogeny.pool") && return :phylogeny_pool
    null_model === :independentswap && return null_model
    null_model === :trialswap && return null_model
    throw(ArgumentError("unsupported phylogenetic null_model: $null_model"))
end

@inline function _tools_sample_labels(n::Integer, labels::Union{Nothing, AbstractVector{<:AbstractString}})
    labels === nothing && return ["sample_$i" for i in 1:n]
    length(labels) == n || throw(ArgumentError("sample_labels length must match number of sample rows"))
    return String.(labels)
end

function _tools_tip_permutation(tree::CompactTree, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}})
    tip_labels === nothing && return collect(1:tree.ntips)
    length(tip_labels) == tree.ntips || throw(ArgumentError("tip_labels length must match number of sample columns"))
    perm = Vector{Int}(undef, tree.ntips)
    for (j, label) in enumerate(tip_labels)
        node = get(tree.tipname_to_id, String(label), Int32(0))
        node == 0 && throw(ArgumentError("sample column label $label is not present in tree tips"))
        idx = findfirst(==(node), tree.tip_ids)
        idx === nothing && throw(ArgumentError("tree tip index not found for label $label"))
        perm[j] = idx
    end
    return perm
end

function _tools_tip_permutation(tree::CompactTree, samp::AbstractDataFrame, tip_labels::Nothing)
    labels = String.(names(samp))
    length(labels) == tree.ntips || throw(ArgumentError("dataframe must have one column per tree tip"))
    return _tools_tip_permutation(tree, labels)
end

function _tools_tip_permutation(tree::CompactTree, samp::AbstractDataFrame, tip_labels::AbstractVector{<:AbstractString})
    length(tip_labels) == ncol(samp) || throw(ArgumentError("tip_labels length must match number of sample columns"))
    return _tools_tip_permutation(tree, tip_labels)
end

function _tools_abundance_only(samp::AbstractMatrix{<:Real}, perm::Vector{Int})
    abund = Matrix{Float64}(undef, size(samp, 1), length(perm))
    @inbounds for i in axes(samp, 1), j in axes(samp, 2)
        abund[i, perm[j]] = Float64(samp[i, j])
    end
    return abund
end

function _tools_presence_only(samp::AbstractMatrix{<:Real}, perm::Vector{Int})
    presence = Matrix{Bool}(undef, size(samp, 1), length(perm))
    @inbounds for i in axes(samp, 1), j in axes(samp, 2)
        presence[i, perm[j]] = samp[i, j] > 0
    end
    return presence
end

function _tools_presence_matrix(samp::AbstractMatrix{<:Real}, perm::Vector{Int})
    presence = Matrix{Bool}(undef, size(samp, 1), length(perm))
    @inbounds for i in axes(samp, 1), j in axes(samp, 2)
        presence[i, perm[j]] = samp[i, j] > 0
    end
    return presence
end

function _tools_tip_paths(tree::CompactTree)
    edge_of_child = zeros(Int32, tree.nnodes)
    @inbounds for edge in 1:tree.nedges
        edge_of_child[Int(tree.child_of_edge[edge])] = Int32(edge)
    end
    paths = Vector{Vector{Int32}}(undef, tree.ntips)
    @inbounds for (i, tip) in enumerate(tree.tip_ids)
        node = Int(tip)
        path = Int32[]
        while node != Int(tree.root)
            edge = edge_of_child[node]
            edge == 0 && break
            push!(path, edge)
            node = Int(tree.parent_of_node[node])
        end
        paths[i] = path
    end
    return paths
end

function _tools_descendant_tip_counts(tree::CompactTree)
    counts = zeros(Int32, tree.nnodes)
    @inbounds for node in tree.postorder
        inode = Int(node)
        if tree.is_tip[inode]
            counts[inode] = 1
        else
            total = Int32(0)
            for child in tree.children[inode]
                total += counts[Int(child)]
            end
            counts[inode] = total
        end
    end
    return counts
end

@inline function _tools_mrca_depth_from_paths(paths::Vector{Vector{Int32}}, present::Vector{Int}, tree::CompactTree)
    length(present) <= 1 && return 0.0
    ref = paths[present[1]]
    shared = 0
    max_shared = length(ref)
    for pos in 1:max_shared
        edge = ref[end - pos + 1]
        ok = true
        for k in 2:length(present)
            path = paths[present[k]]
            if pos > length(path) || path[end - pos + 1] != edge
                ok = false
                break
            end
        end
        ok || break
        shared = pos
    end
    shared == 0 && return 0.0
    mrca_edge = ref[end - shared + 1]
    return tree.dist_from_root[Int(tree.child_of_edge[Int(mrca_edge)])]
end

@inline function _tools_update_rank!(rank_count::Vector{Float64}, obs::Vector{Float64}, rand::Vector{Float64})
    @inbounds for i in eachindex(obs)
        if isfinite(obs[i]) && isfinite(rand[i])
            if rand[i] < obs[i]
                rank_count[i] += 1.0
            elseif rand[i] == obs[i]
                rank_count[i] += 0.5
            end
        end
    end
    return rank_count
end

function _randomize_richness!(dest::Matrix{Bool}, richness::Vector{Int}, rng::AbstractRNG)
    fill!(dest, false)
    nspecies = size(dest, 2)
    @inbounds for i in axes(dest, 1)
        r = richness[i]
        r <= 0 && continue
        cols = randperm(rng, nspecies)
        for k in 1:r
            dest[i, cols[k]] = true
        end
    end
    return dest
end

function _randomize_frequency!(dest::Matrix{Bool}, freq::Vector{Int}, rng::AbstractRNG)
    fill!(dest, false)
    nsamples = size(dest, 1)
    @inbounds for j in axes(dest, 2)
        f = freq[j]
        f <= 0 && continue
        rows = randperm(rng, nsamples)
        for k in 1:f
            dest[rows[k], j] = true
        end
    end
    return dest
end

function _randomize_taxa_labels!(dest::Matrix{Bool}, presence::Matrix{Bool}, rng::AbstractRNG)
    perm = randperm(rng, size(presence, 2))
    @inbounds for i in axes(presence, 1), j in axes(presence, 2)
        dest[i, j] = presence[i, perm[j]]
    end
    return dest
end

@inline function _swap_checkerboard!(m::Matrix{Bool}, i::Int, j::Int, k::Int, l::Int)
    a = m[i, k]
    b = m[i, l]
    c = m[j, k]
    d = m[j, l]
    if (a && d && !b && !c) || (!a && !d && b && c)
        m[i, k] = c
        m[j, k] = a
        m[i, l] = d
        m[j, l] = b
        return true
    end
    return false
end

function _random_pair(rng::AbstractRNG, n::Int)
    n >= 2 || throw(ArgumentError("swap null models require at least two rows and two columns"))
    a = rand(rng, 1:n)
    b = rand(rng, 1:(n - 1))
    b >= a && (b += 1)
    return a, b
end

function _randomize_trialswap!(dest::Matrix{Bool}, presence::Matrix{Bool}, iterations::Integer, rng::AbstractRNG)
    copyto!(dest, presence)
    nrow, ncol = size(dest)
    (nrow >= 2 && ncol >= 2) || return dest
    @inbounds for _ in 1:iterations
        i, j = _random_pair(rng, nrow)
        k, l = _random_pair(rng, ncol)
        _swap_checkerboard!(dest, i, j, k, l)
    end
    return dest
end

function _randomize_independentswap!(
    dest::Matrix{Bool},
    presence::Matrix{Bool},
    iterations::Integer,
    rng::AbstractRNG;
    max_attempts_per_swap::Integer = 10_000,
)
    copyto!(dest, presence)
    nrow, ncol = size(dest)
    (nrow >= 2 && ncol >= 2) || return dest
    @inbounds for _ in 1:iterations
        swapped = false
        for _attempt in 1:max_attempts_per_swap
            i, j = _random_pair(rng, nrow)
            k, l = _random_pair(rng, ncol)
            if _swap_checkerboard!(dest, i, j, k, l)
                swapped = true
                break
            end
        end
        swapped || break
    end
    return dest
end

function _tools_mpd_site!(
    tree::CompactTree,
    selected::AbstractVector{Bool},
    weights::Union{Nothing, AbstractVector{Float64}},
    node_mass::Vector{Float64},
)
    fill!(node_mass, 0.0)
    if weights === nothing
        @inbounds for (i, tip) in enumerate(tree.tip_ids)
            node_mass[Int(tip)] = selected[i] ? 1.0 : 0.0
        end
    else
        @inbounds for (i, tip) in enumerate(tree.tip_ids)
            node_mass[Int(tip)] = selected[i] ? weights[i] : 0.0
        end
    end
    @inbounds for node in tree.postorder
        inode = Int(node)
        tree.is_tip[inode] && continue
        total = 0.0
        for child in tree.children[inode]
            total += node_mass[Int(child)]
        end
        node_mass[inode] = total
    end
    total_mass = node_mass[Int(tree.root)]
    total_mass <= 0.0 && return NaN
    if weights === nothing
        total_mass < 2.0 && return NaN
        dist_sum = 0.0
        @inbounds for edge in 1:tree.nedges
            child = Int(tree.child_of_edge[edge])
            left = node_mass[child]
            dist_sum += tree.edge_length[edge] * left * (total_mass - left)
        end
        return 2.0 * dist_sum / (total_mass * (total_mass - 1.0))
    else
        dist_sum = 0.0
        @inbounds for edge in 1:tree.nedges
            child = Int(tree.child_of_edge[edge])
            left = node_mass[child]
            dist_sum += tree.edge_length[edge] * left * (total_mass - left)
        end
        return 2.0 * dist_sum / (total_mass * total_mass)
    end
end
