Base.@kwdef struct PDResult
    sample_labels::Vector{String} = String[]
    pd::Vector{Float64} = Float64[]
    richness::Vector{Int} = Int[]
    include_root::Bool = true
end

Base.@kwdef struct SESPDResult
    sample_labels::Vector{String} = String[]
    ntaxa::Vector{Int} = Int[]
    pd_obs::Vector{Float64} = Float64[]
    pd_rand_mean::Vector{Float64} = Float64[]
    pd_rand_sd::Vector{Float64} = Float64[]
    pd_obs_rank::Vector{Float64} = Float64[]
    pd_obs_z::Vector{Float64} = Float64[]
    pd_obs_p::Vector{Float64} = Float64[]
    runs::Int = 0
    null_model::Symbol = :unknown
    include_root::Bool = true
end

function _pd_from_presence!(out::Vector{Float64}, richness::Vector{Int}, presence::AbstractMatrix{Bool}, tree::CompactTree, paths::Vector{Vector{Int32}}, stamp::Vector{Int32}, include_root::Bool)
    fill!(out, 0.0)
    fill!(richness, 0)
    present = Int[]
    @inbounds for i in axes(presence, 1)
        empty!(present)
        for j in axes(presence, 2)
            presence[i, j] && push!(present, j)
        end
        richness[i] = length(present)
        if isempty(present)
            out[i] = 0.0
            continue
        elseif length(present) == 1 && !include_root
            out[i] = NaN
            continue
        end
        marker = Int32(i)
        total = 0.0
        for tip_idx in present
            for edge32 in paths[tip_idx]
                edge = Int(edge32)
                if stamp[edge] != marker
                    stamp[edge] = marker
                    total += tree.edge_length[edge]
                end
            end
        end
        !include_root && (total -= _tools_mrca_depth_from_paths(paths, present, tree))
        out[i] = total
    end
    return out
end

function pd(samp::AbstractMatrix{<:Real}, tree::CompactTree; include_root::Bool = true, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing)
    size(samp, 2) == tree.ntips || throw(ArgumentError("samp must have one column per tree tip"))
    perm = _tools_tip_permutation(tree, tip_labels)
    presence = _tools_presence_matrix(samp, perm)
    paths = _tools_tip_paths(tree)
    stamp = zeros(Int32, tree.nedges)
    pd_values = Vector{Float64}(undef, size(samp, 1))
    richness = Vector{Int}(undef, size(samp, 1))
    _pd_from_presence!(pd_values, richness, presence, tree, paths, stamp, include_root)
    return PDResult(sample_labels = _tools_sample_labels(size(samp, 1), sample_labels), pd = pd_values, richness = richness, include_root = include_root)
end

function pd_dataframe(samp::AbstractMatrix{<:Real}, tree::CompactTree; kwargs...)
    result = pd(samp, tree; kwargs...)
    return DataFrame(sample = result.sample_labels, PD = result.pd, SR = result.richness)
end

function ses_pd(samp::AbstractMatrix{<:Real}, tree::CompactTree; null_model::Symbol = :taxa_labels, runs::Integer = 999, iterations::Integer = 1000, include_root::Bool = true, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, rng::AbstractRNG = Random.default_rng(), max_attempts_per_swap::Integer = 10_000)
    runs >= 1 || throw(ArgumentError("runs must be positive"))
    model = _tools_normalize_phylo_null(null_model)
    size(samp, 2) == tree.ntips || throw(ArgumentError("samp must have one column per tree tip"))
    perm = _tools_tip_permutation(tree, tip_labels)
    presence = _tools_presence_matrix(samp, perm)

    paths = _tools_tip_paths(tree)
    stamp = zeros(Int32, tree.nedges)
    nsites = size(presence, 1)
    obs_pd = Vector{Float64}(undef, nsites)
    richness = Vector{Int}(undef, nsites)
    _pd_from_presence!(obs_pd, richness, presence, tree, paths, stamp, include_root)

    freq = Int.(vec(sum(presence; dims = 1)))
    rand_presence = similar(presence)
    tmp_presence = similar(presence)
    rand_pd = Vector{Float64}(undef, nsites)
    rand_richness = Vector{Int}(undef, nsites)
    sum_pd = zeros(Float64, nsites)
    sumsq_pd = zeros(Float64, nsites)
    rank_count = ones(Float64, nsites)

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
            copyto!(tmp_presence, rand_presence)
            _randomize_taxa_labels!(rand_presence, tmp_presence, rng)
        end
        _pd_from_presence!(rand_pd, rand_richness, rand_presence, tree, paths, stamp, include_root)
        @inbounds for i in 1:nsites
            if isfinite(rand_pd[i])
                sum_pd[i] += rand_pd[i]
                sumsq_pd[i] += rand_pd[i] * rand_pd[i]
            end
        end
        _tools_update_rank!(rank_count, obs_pd, rand_pd)
    end

    means = sum_pd ./ Float64(runs)
    sds = Vector{Float64}(undef, nsites)
    z = Vector{Float64}(undef, nsites)
    pvals = Vector{Float64}(undef, nsites)
    ranks = rank_count
    @inbounds for i in 1:nsites
        var = runs > 1 ? (sumsq_pd[i] - runs * means[i] * means[i]) / (runs - 1) : NaN
        sds[i] = var >= 0.0 ? sqrt(max(var, 0.0)) : NaN
        z[i] = isfinite(sds[i]) && sds[i] > 0.0 ? (obs_pd[i] - means[i]) / sds[i] : NaN
        pvals[i] = ranks[i] / (runs + 1)
    end

    return SESPDResult(sample_labels = _tools_sample_labels(nsites, sample_labels), ntaxa = richness, pd_obs = obs_pd, pd_rand_mean = means, pd_rand_sd = sds, pd_obs_rank = ranks, pd_obs_z = z, pd_obs_p = pvals, runs = Int(runs), null_model = model, include_root = include_root)
end

function ses_pd_dataframe(samp::AbstractMatrix{<:Real}, tree::CompactTree; kwargs...)
    result = ses_pd(samp, tree; kwargs...)
    return DataFrame(sample = result.sample_labels, ntaxa = result.ntaxa, pd_obs = result.pd_obs, pd_rand_mean = result.pd_rand_mean, pd_rand_sd = result.pd_rand_sd, pd_obs_rank = result.pd_obs_rank, pd_obs_z = result.pd_obs_z, pd_obs_p = result.pd_obs_p, runs = fill(result.runs, length(result.pd_obs)))
end
