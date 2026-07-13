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
