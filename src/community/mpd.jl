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

"""
    ses_mpd(community, tree; null_model=:taxa_labels, runs=100,
            abundance_weighted=false, kwargs...)

MPD standardized effect size using sparse null communities and reusable virtual
tree workspaces. Abundance weighting retains the original tip-position weights
while randomizing presence. Null tips with zero original weight have zero mass.
Finite draws determine online mean and sample SD; no finite draws or zero SD
produce NaN z. Null-model parameters and rank p-values follow ses_pd.
"""
function ses_mpd(samp::AbstractMatrix{<:Real}, tree::CompactTree; null_model::Symbol = :taxa_labels, runs::Integer = 100, iterations::Integer = 1000, abundance_weighted::Bool = false, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing, rng::AbstractRNG = Random.default_rng(), max_attempts_per_swap::Integer = 10_000)
    result = _ses_summary(:mpd, samp, tree; null_model, runs, iterations, abundance_weighted,
        tip_labels, sample_labels, rng, max_attempts_per_swap)
    return SESMPDResult(sample_labels = result.labels, ntaxa = result.richness,
        mpd_obs = result.obs, mpd_rand_mean = result.means, mpd_rand_sd = result.sds,
        mpd_obs_rank = result.ranks, mpd_obs_z = result.z, mpd_nri = -result.z,
        mpd_obs_p = result.p, runs = Int(runs), null_model = result.model,
        abundance_weighted = abundance_weighted)
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
