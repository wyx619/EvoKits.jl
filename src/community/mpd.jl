Base.@kwdef struct SESMPDResult
    sample_labels::Vector{String} = String[]
    ntaxa::Vector{Int} = Int[]
    mpd_obs::Vector{Float64} = Float64[]
    mpd_rand_mean::Vector{Float64} = Float64[]
    mpd_rand_sd::Vector{Float64} = Float64[]
    mpd_obs_rank::Vector{Float64} = Float64[]
    mpd_obs_z::Vector{Float64} = Float64[]
    mpd_nri::Vector{Float64} = Float64[]
    mpd_obs_p::Vector{Float64} = Float64[]
    runs::Int = 0
    null_model::Symbol = :unknown
    abundance_weighted::Bool = false
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

function ses_mpd(samp::AbstractMatrix{<:Real}, tree::CompactTree; null_model::Symbol = :taxa_labels, runs::Integer = 100, iterations::Integer = 1000, abundance_weighted::Bool = false, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, rng::AbstractRNG = Random.default_rng(), max_attempts_per_swap::Integer = 10_000)
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
    node_mass = zeros(Float64, tree.nnodes)
    @inbounds for i in 1:nsites
        if abundance_weighted
            obs[i] = _tools_mpd_site!(tree, @view(presence[i, :]), @view(abund[i, :]), node_mass)
        else
            selected = findall(@view(presence[i, :]))
            obs[i] = length(selected) >= 2 ? first(_community_metrics(tree, selected, query_index)) : NaN
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
                rand_obs[i] = _tools_mpd_site!(tree, @view(rand_presence[i, :]), @view(abund[i, :]), node_mass)
            else
                selected = findall(@view(rand_presence[i, :]))
                rand_obs[i] = length(selected) >= 2 ? first(_community_metrics(tree, selected, query_index)) : NaN
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
    nri = Vector{Float64}(undef, nsites)
    pv = Vector{Float64}(undef, nsites)
    ranks = rank_count
    @inbounds for i in 1:nsites
        var = moments_n[i] > 1 ? moments_m2[i] / (moments_n[i] - 1) : NaN
        sdv[i] = var >= 0.0 ? sqrt(max(var, 0.0)) : NaN
        zv[i] = isfinite(sdv[i]) && sdv[i] > 0.0 ? (obs[i] - meanv[i]) / sdv[i] : NaN
        nri[i] = isfinite(zv[i]) ? -zv[i] : NaN
        pv[i] = ranks[i] / (runs + 1)
    end
    return SESMPDResult(sample_labels = _tools_sample_labels(nsites, sample_labels), ntaxa = collect(richness), mpd_obs = obs, mpd_rand_mean = meanv, mpd_rand_sd = sdv, mpd_obs_rank = ranks, mpd_obs_z = zv, mpd_nri = nri, mpd_obs_p = pv, runs = Int(runs), null_model = model, abundance_weighted = abundance_weighted)
end

function ses_mpd(samp::AbstractDataFrame, tree::CompactTree; kwargs...)
    labels = haskey(kwargs, :tip_labels) ? kwargs[:tip_labels] : String.(names(samp))
    return ses_mpd(Matrix{Float64}(samp), tree; kwargs..., tip_labels = labels)
end

function ses_mpd_dataframe(samp::AbstractMatrix{<:Real}, tree::CompactTree; kwargs...)
    result = ses_mpd(samp, tree; kwargs...)
    return DataFrame(sample = result.sample_labels, ntaxa = result.ntaxa, mpd_obs = result.mpd_obs, mpd_rand_mean = result.mpd_rand_mean, mpd_rand_sd = result.mpd_rand_sd, mpd_obs_rank = result.mpd_obs_rank, mpd_obs_z = result.mpd_obs_z, mpd_nri = result.mpd_nri, mpd_obs_p = result.mpd_obs_p, runs = fill(result.runs, length(result.mpd_obs)))
end

function ses_mpd_dataframe(samp::AbstractDataFrame, tree::CompactTree; kwargs...)
    result = ses_mpd(samp, tree; kwargs...)
    return DataFrame(sample = result.sample_labels, ntaxa = result.ntaxa, mpd_obs = result.mpd_obs, mpd_rand_mean = result.mpd_rand_mean, mpd_rand_sd = result.mpd_rand_sd, mpd_obs_rank = result.mpd_obs_rank, mpd_obs_z = result.mpd_obs_z, mpd_nri = result.mpd_nri, mpd_obs_p = result.mpd_obs_p, runs = fill(result.runs, length(result.mpd_obs)))
end
