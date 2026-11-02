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

"""
    ses_mntd(community, tree; null_model=:taxa_labels, runs=100,
             abundance_weighted=false, kwargs...)

MNTD standardized effect size using sparse null communities and reusable virtual
tree workspaces. Abundance weights remain at original tree-tip positions;
randomized tips with zero weight remain nearest-neighbor candidates. Finite
draws determine online mean and sample SD; no finite draws or zero SD produce
NaN z. Null-model parameters and rank p-values follow ses_pd.
"""
function ses_mntd(samp::AbstractMatrix{<:Real}, tree::CompactTree; null_model::Symbol = :taxa_labels, runs::Integer = 100, iterations::Integer = 1000, abundance_weighted::Bool = false, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, rng::AbstractRNG = Random.default_rng(), max_attempts_per_swap::Integer = 10_000)
    result = _ses_summary(:mntd, samp, tree; null_model, runs, iterations, abundance_weighted,
        tip_labels, sample_labels, rng, max_attempts_per_swap)
    return SESMNTDResult(sample_labels = result.labels, ntaxa = result.richness,
        mntd_obs = result.obs, mntd_rand_mean = result.means, mntd_rand_sd = result.sds,
        mntd_obs_rank = result.ranks, mntd_obs_z = result.z, mntd_nti = -result.z,
        mntd_obs_p = result.p, runs = Int(runs), null_model = result.model,
        abundance_weighted = abundance_weighted)
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
