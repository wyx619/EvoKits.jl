Base.@kwdef struct SESMNTDResult
    sample_labels::Vector{String} = String[]
    ntaxa::Vector{Int} = Int[]
    mntd_obs::Vector{Float64} = Float64[]
    mntd_rand_mean::Vector{Float64} = Float64[]
    mntd_rand_sd::Vector{Float64} = Float64[]
    mntd_obs_rank::Vector{Float64} = Float64[]
    mntd_obs_z::Vector{Float64} = Float64[]
    mntd_nti::Vector{Float64} = Float64[]
    mntd_obs_p::Vector{Float64} = Float64[]
    runs::Int = 0
    null_model::Symbol = :unknown
    abundance_weighted::Bool = false
end

@inline function _mntd_merge_top2!(d1::Float64, t1::Int32, d2::Float64, t2::Int32, nd1::Float64, nt1::Int32, nd2::Float64, nt2::Int32)
    if nt1 != 0
        if nt1 == t1
            nd1 < d1 && (d1 = nd1)
        elseif nt1 == t2
            nd1 < d2 && (d2 = nd1)
        elseif nd1 < d1
            d2, t2 = d1, t1
            d1, t1 = nd1, nt1
        elseif nd1 < d2
            d2, t2 = nd1, nt1
        end
    end
    if nt2 != 0
        if nt2 == t1
            nd2 < d1 && (d1 = nd2)
        elseif nt2 == t2
            nd2 < d2 && (d2 = nd2)
        elseif nd2 < d1
            d2, t2 = d1, t1
            d1, t1 = nd2, nt2
        elseif nd2 < d2
            d2, t2 = nd2, nt2
        end
    end
    if d2 < d1
        d1, d2 = d2, d1
        t1, t2 = t2, t1
    end
    return d1, t1, d2, t2
end

@inline function _mntd_edge_length(tree::CompactTree, node::Int)
    parent = Int(tree.parent_of_node[node])
    parent == 0 && return 0.0
    return tree.dist_from_root[node] - tree.dist_from_root[parent]
end

function _mntd_site!(tree::CompactTree, selected::AbstractVector{Bool}, tip_index_of_node::Vector{Int32}, weights::Union{Nothing, AbstractVector{Float64}}, down1d::Vector{Float64}, down1t::Vector{Int32}, down2d::Vector{Float64}, down2t::Vector{Int32}, up1d::Vector{Float64}, up1t::Vector{Int32}, up2d::Vector{Float64}, up2t::Vector{Int32}, pref1d::Vector{Float64}, pref1t::Vector{Int32}, pref2d::Vector{Float64}, pref2t::Vector{Int32}, suff1d::Vector{Float64}, suff1t::Vector{Int32}, suff2d::Vector{Float64}, suff2t::Vector{Int32}, child1d::Vector{Float64}, child1t::Vector{Int32}, child2d::Vector{Float64}, child2t::Vector{Int32})
    fill!(down1d, Inf); fill!(down1t, Int32(0)); fill!(down2d, Inf); fill!(down2t, Int32(0))
    fill!(up1d, Inf); fill!(up1t, Int32(0)); fill!(up2d, Inf); fill!(up2t, Int32(0))
    @inbounds for node in tree.postorder
        inode = Int(node)
        if tree.is_tip[inode]
            tip_i = Int(tip_index_of_node[inode])
            if tip_i != 0 && selected[tip_i]
                down1d[inode] = 0.0
                down1t[inode] = Int32(inode)
            end
            continue
        end
        b1d = Inf; b1t = Int32(0); b2d = Inf; b2t = Int32(0)
        for child in tree.children[inode]
            c = Int(child)
            len = _mntd_edge_length(tree, c)
            b1d, b1t, b2d, b2t = _mntd_merge_top2!(b1d, b1t, b2d, b2t, down1d[c] + len, down1t[c], down2d[c] + len, down2t[c])
        end
        down1d[inode] = b1d; down1t[inode] = b1t; down2d[inode] = b2d; down2t[inode] = b2t
    end

    @inbounds for node in tree.preorder
        inode = Int(node)
        children = tree.children[inode]
        d = length(children)
        d == 0 && continue
        @inbounds for j in 1:d
            c = Int(children[j])
            child1d[j] = down1d[c] + _mntd_edge_length(tree, c)
            child1t[j] = down1t[c]
            child2d[j] = down2d[c] + _mntd_edge_length(tree, c)
            child2t[j] = down2t[c]
        end
        pref1d[1] = up1d[inode]; pref1t[1] = up1t[inode]; pref2d[1] = up2d[inode]; pref2t[1] = up2t[inode]
        @inbounds for j in 1:d
            pref1d[j + 1], pref1t[j + 1], pref2d[j + 1], pref2t[j + 1] = _mntd_merge_top2!(pref1d[j], pref1t[j], pref2d[j], pref2t[j], child1d[j], child1t[j], child2d[j], child2t[j])
        end
        suff1d[d + 1] = Inf; suff1t[d + 1] = Int32(0); suff2d[d + 1] = Inf; suff2t[d + 1] = Int32(0)
        @inbounds for j in d:-1:1
            suff1d[j], suff1t[j], suff2d[j], suff2t[j] = _mntd_merge_top2!(child1d[j], child1t[j], child2d[j], child2t[j], suff1d[j + 1], suff1t[j + 1], suff2d[j + 1], suff2t[j + 1])
        end
        @inbounds for j in 1:d
            c = Int(children[j])
            best1d, best1t, best2d, best2t = _mntd_merge_top2!(pref1d[j], pref1t[j], pref2d[j], pref2t[j], suff1d[j + 1], suff1t[j + 1], suff2d[j + 1], suff2t[j + 1])
            len = _mntd_edge_length(tree, c)
            up1d[c] = best1d + len; up1t[c] = best1t; up2d[c] = best2d + len; up2t[c] = best2t
        end
    end

    total = 0.0
    weight_sum = 0.0
    found = 0
    @inbounds for (i, tip) in enumerate(tree.tip_ids)
        selected[i] || continue
        node = Int(tip)
        d1 = up1d[node]; t1 = up1t[node]; d2 = up2d[node]
        nearest = (t1 != node && isfinite(d1)) ? d1 : d2
        if !isfinite(nearest)
            continue
        end
        if weights === nothing
            total += nearest
            found += 1
        else
            w = weights[i]
            weight_sum += w
            total += w * nearest
        end
    end
    if weights === nothing
        found < 1 && return NaN
        return total / Float64(found)
    end
    weight_sum <= 0.0 && return NaN
    return total / weight_sum
end

function ses_mntd(samp::AbstractMatrix{<:Real}, tree::CompactTree; null_model::Symbol = :taxa_labels, runs::Integer = 100, iterations::Integer = 1000, abundance_weighted::Bool = false, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, rng::AbstractRNG = Random.default_rng(), max_attempts_per_swap::Integer = 10_000)
    runs >= 1 || throw(ArgumentError("runs must be positive"))
    iterations >= 0 || throw(ArgumentError("iterations must be non-negative"))
    model = _tools_normalize_phylo_null(null_model)
    size(samp, 2) == tree.ntips || throw(ArgumentError("samp must have one column per tree tip"))
    perm = _tools_tip_permutation(tree, tip_labels)
    presence = _tools_presence_only(samp, perm)
    abund = abundance_weighted ? _tools_abundance_only(samp, perm) : nothing
    nsites = size(presence, 1)
    query_index = abundance_weighted ? nothing : _build_tree_query_index(tree)
    obs = Vector{Float64}(undef, nsites)
    tip_index_of_node = zeros(Int32, tree.nnodes)
    @inbounds for (i, tip) in enumerate(tree.tip_ids)
        tip_index_of_node[Int(tip)] = Int32(i)
    end
    maxdeg = maximum(length.(tree.children))
    down1d = fill(Inf, tree.nnodes); down1t = fill(Int32(0), tree.nnodes); down2d = fill(Inf, tree.nnodes); down2t = fill(Int32(0), tree.nnodes)
    up1d = fill(Inf, tree.nnodes); up1t = fill(Int32(0), tree.nnodes); up2d = fill(Inf, tree.nnodes); up2t = fill(Int32(0), tree.nnodes)
    pref1d = fill(Inf, maxdeg + 1); pref1t = fill(Int32(0), maxdeg + 1); pref2d = fill(Inf, maxdeg + 1); pref2t = fill(Int32(0), maxdeg + 1)
    suff1d = fill(Inf, maxdeg + 1); suff1t = fill(Int32(0), maxdeg + 1); suff2d = fill(Inf, maxdeg + 1); suff2t = fill(Int32(0), maxdeg + 1)
    child1d = fill(Inf, maxdeg); child1t = fill(Int32(0), maxdeg); child2d = fill(Inf, maxdeg); child2t = fill(Int32(0), maxdeg)
    @inbounds for i in 1:nsites
        if abundance_weighted
            obs[i] = _mntd_site!(tree, @view(presence[i, :]), tip_index_of_node, @view(abund[i, :]), down1d, down1t, down2d, down2t, up1d, up1t, up2d, up2t, pref1d, pref1t, pref2d, pref2t, suff1d, suff1t, suff2d, suff2t, child1d, child1t, child2d, child2t)
        else
            selected = findall(@view(presence[i, :]))
            obs[i] = length(selected) >= 2 ? last(_community_metrics(tree, selected, query_index)) : NaN
        end
    end
    richness = vec(sum(presence; dims = 2))
    rand_presence = similar(presence)
    rand_obs = Vector{Float64}(undef, nsites)
    moments_n = zeros(Int, nsites)
    meanv = zeros(Float64, nsites)
    moments_m2 = zeros(Float64, nsites)
    rank_count = ones(Float64, nsites)
    freq = Int.(vec(sum(presence; dims = 1)))
    for _ in 1:runs
        if model === :taxa_labels
            _randomize_taxa_labels!(rand_presence, presence, rng)
        elseif model === :richness || model === :sample_pool || model === :phylogeny_pool
            _randomize_richness!(rand_presence, richness, rng)
        elseif model === :frequency
            _randomize_frequency!(rand_presence, freq, rng)
        elseif model === :trialswap
            _randomize_trialswap!(rand_presence, presence, iterations, rng)
        elseif model === :independentswap
            _randomize_independentswap!(rand_presence, presence, iterations, rng; max_attempts_per_swap = max_attempts_per_swap)
        end
        if model === :phylogeny_pool
            tmp = similar(rand_presence)
            copyto!(tmp, rand_presence)
            _randomize_taxa_labels!(rand_presence, tmp, rng)
        end
        @inbounds for i in 1:nsites
            if abundance_weighted
                rand_obs[i] = _mntd_site!(tree, @view(rand_presence[i, :]), tip_index_of_node, @view(abund[i, :]), down1d, down1t, down2d, down2t, up1d, up1t, up2d, up2t, pref1d, pref1t, pref2d, pref2t, suff1d, suff1t, suff2d, suff2t, child1d, child1t, child2d, child2t)
            else
                selected = findall(@view(rand_presence[i, :]))
                rand_obs[i] = length(selected) >= 2 ? last(_community_metrics(tree, selected, query_index)) : NaN
            end
        end
        @inbounds for i in 1:nsites
            if isfinite(rand_obs[i])
                moments_n[i] += 1
                delta = rand_obs[i] - meanv[i]
                meanv[i] += delta / moments_n[i]
                moments_m2[i] += delta * (rand_obs[i] - meanv[i])
            end
        end
        _tools_update_rank!(rank_count, obs, rand_obs)
    end
    sdv = Vector{Float64}(undef, nsites)
    zv = Vector{Float64}(undef, nsites)
    nti = Vector{Float64}(undef, nsites)
    pv = Vector{Float64}(undef, nsites)
    ranks = rank_count
    @inbounds for i in 1:nsites
        var = moments_n[i] > 1 ? moments_m2[i] / (moments_n[i] - 1) : NaN
        sdv[i] = var >= 0.0 ? sqrt(max(var, 0.0)) : NaN
        zv[i] = isfinite(sdv[i]) && sdv[i] > 0.0 ? (obs[i] - meanv[i]) / sdv[i] : NaN
        nti[i] = isfinite(zv[i]) ? -zv[i] : NaN
        pv[i] = ranks[i] / (runs + 1)
    end
    return SESMNTDResult(sample_labels = _tools_sample_labels(nsites, sample_labels), ntaxa = collect(richness), mntd_obs = obs, mntd_rand_mean = meanv, mntd_rand_sd = sdv, mntd_obs_rank = ranks, mntd_obs_z = zv, mntd_nti = nti, mntd_obs_p = pv, runs = Int(runs), null_model = model, abundance_weighted = abundance_weighted)
end

function ses_mntd(samp::AbstractDataFrame, tree::CompactTree; kwargs...)
    labels = haskey(kwargs, :tip_labels) ? kwargs[:tip_labels] : String.(names(samp))
    return ses_mntd(Matrix{Float64}(samp), tree; kwargs..., tip_labels = labels)
end

function ses_mntd_dataframe(samp::AbstractMatrix{<:Real}, tree::CompactTree; kwargs...)
    result = ses_mntd(samp, tree; kwargs...)
    return DataFrame(sample = result.sample_labels, ntaxa = result.ntaxa, mntd_obs = result.mntd_obs, mntd_rand_mean = result.mntd_rand_mean, mntd_rand_sd = result.mntd_rand_sd, mntd_obs_rank = result.mntd_obs_rank, mntd_obs_z = result.mntd_obs_z, mntd_nti = result.mntd_nti, mntd_obs_p = result.mntd_obs_p, runs = fill(result.runs, length(result.mntd_obs)))
end

function ses_mntd_dataframe(samp::AbstractDataFrame, tree::CompactTree; kwargs...)
    result = ses_mntd(samp, tree; kwargs...)
    return DataFrame(sample = result.sample_labels, ntaxa = result.ntaxa, mntd_obs = result.mntd_obs, mntd_rand_mean = result.mntd_rand_mean, mntd_rand_sd = result.mntd_rand_sd, mntd_obs_rank = result.mntd_obs_rank, mntd_obs_z = result.mntd_obs_z, mntd_nti = result.mntd_nti, mntd_obs_p = result.mntd_obs_p, runs = fill(result.runs, length(result.mntd_obs)))
end
