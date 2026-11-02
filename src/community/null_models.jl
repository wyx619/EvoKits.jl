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

function _random_pair(rng::AbstractRNG, n::Int)
    n >= 2 || throw(ArgumentError("swap null models require at least two rows and two columns"))
    a = rand(rng, 1:n)
    b = rand(rng, 1:(n - 1))
    b >= a && (b += 1)
    return a, b
end

"""Reusable sorted presence rows and O(nsites + nspecies) permutation buffers."""
struct _CommunityNullWorkspace
    rows::Vector{Vector{Int32}}
    richness::Vector{Int}
    frequency::Vector{Int}
    species_permutation::Vector{Int}
    inverse_permutation::Vector{Int}
    site_permutation::Vector{Int}
end

function _CommunityNullWorkspace(observed, nspecies)
    frequency = zeros(Int, nspecies)
    for row in observed, tip in row
        frequency[tip] += 1
    end
    rows = [sizehint!(Int32[], length(row)) for row in observed]
    return _CommunityNullWorkspace(rows, length.(observed), frequency,
        Vector{Int}(undef, nspecies), Vector{Int}(undef, nspecies),
        Vector{Int}(undef, length(observed)))
end

function _randomize_richness!(workspace::_CommunityNullWorkspace, rng)
    for (site, row) in enumerate(workspace.rows)
        richness = workspace.richness[site]
        resize!(row, richness)
        richness == 0 && continue
        # Reuse the full shuffle buffer to retain historical RNG consumption.
        # Storage is linear in species, never sites times species.
        randperm!(rng, workspace.species_permutation)
        for k in 1:richness
            row[k] = Int32(workspace.species_permutation[k])
        end
        sort!(row)
    end
    return workspace.rows
end

function _randomize_frequency!(workspace::_CommunityNullWorkspace, rng)
    foreach(empty!, workspace.rows)
    for (tip, frequency) in enumerate(workspace.frequency)
        frequency == 0 && continue
        randperm!(rng, workspace.site_permutation)
        for k in 1:frequency
            push!(workspace.rows[workspace.site_permutation[k]], Int32(tip))
        end
    end
    return workspace.rows
end

function _randomize_taxa_labels!(workspace::_CommunityNullWorkspace, source, rng)
    randperm!(rng, workspace.species_permutation)
    for (target, origin) in enumerate(workspace.species_permutation)
        workspace.inverse_permutation[origin] = target
    end
    # source may be workspace.rows for phylogeny_pool. Each row is remapped
    # in place after the common permutation has been constructed.
    for site in eachindex(source)
        row, observed = workspace.rows[site], source[site]
        resize!(row, length(observed))
        for k in eachindex(observed)
            row[k] = Int32(workspace.inverse_permutation[observed[k]])
        end
        sort!(row)
    end
    return workspace.rows
end

@inline function _ses_contains(row, tip)
    index = searchsortedfirst(row, tip)
    return index <= length(row) && row[index] == tip
end

function _ses_replace_tip!(row, old, new)
    deleteat!(row, searchsortedfirst(row, old))
    insert!(row, searchsortedfirst(row, new), Int32(new))
    return row
end

function _swap_checkerboard!(rows::Vector{Vector{Int32}}, i, j, k, l)
    left, right = rows[i], rows[j]
    a, b = _ses_contains(left, k), _ses_contains(left, l)
    c, d = _ses_contains(right, k), _ses_contains(right, l)
    if a && d && !b && !c
        _ses_replace_tip!(left, k, l)
        _ses_replace_tip!(right, l, k)
    elseif !a && !d && b && c
        _ses_replace_tip!(left, l, k)
        _ses_replace_tip!(right, k, l)
    else
        return false
    end
    return true
end

function _randomize_swaps!(workspace, observed, model, iterations, max_attempts, rng)
    for (row, source) in zip(workspace.rows, observed)
        resize!(row, length(source))
        copyto!(row, source)
    end
    nrow, ncol = length(observed), length(workspace.frequency)
    (nrow >= 2 && ncol >= 2) || return workspace.rows
    attempts = model === :trialswap ? 1 : max_attempts
    for _ in 1:iterations
        swapped = false
        for _ in 1:attempts
            i, j = _random_pair(rng, nrow)
            k, l = _random_pair(rng, ncol)
            if _swap_checkerboard!(workspace.rows, i, j, k, l)
                swapped = true
                break
            end
        end
        !swapped && model === :independentswap && break
    end
    return workspace.rows
end

function _randomize_community!(workspace, observed, model, iterations, max_attempts, rng)
    if model === :taxa_labels
        return _randomize_taxa_labels!(workspace, observed, rng)
    elseif model === :frequency
        return _randomize_frequency!(workspace, rng)
    elseif model in (:richness, :sample_pool, :phylogeny_pool)
        _randomize_richness!(workspace, rng)
        model === :phylogeny_pool && _randomize_taxa_labels!(workspace, workspace.rows, rng)
        return workspace.rows
    end
    return _randomize_swaps!(workspace, observed, model, iterations, max_attempts, rng)
end
