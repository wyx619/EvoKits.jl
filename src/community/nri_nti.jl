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

function _nri_tip_positions(tree::CompactTree, values, name::AbstractString)
    output = if eltype(values) <: Integer
        Int.(values)
    else
        positions = Int[]
        for label0 in values
            label = String(label0)
            node = get(tree.tipname_to_id, label, Int32(0))
            node == 0 && throw(ArgumentError("$name contains unknown tree tip: $label"))
            position = findfirst(==(node), tree.tip_ids)
            position === nothing && throw(ArgumentError("$name tip is not indexed: $label"))
            push!(positions, position)
        end
        positions
    end
    all((1 .<= output) .& (output .<= tree.ntips)) || throw(ArgumentError("$name contains tip positions outside 1:$(tree.ntips)"))
    length(unique(output)) == length(output) || throw(ArgumentError("$name contains duplicate tips"))
    return output
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

"""
    nri_nti(tree, pool, observed; nreps = 999, seed = 20260713)

Compute NRI and NTI for an observed community sampled without replacement from
its local `pool`. `pool` and `observed` may contain tip positions or tip labels.
The observed community must contain at least two tips and be a subset of the
pool. Results are deterministic for the same `seed`.
"""
function nri_nti(
    tree::CompactTree,
    pool,
    observed;
    nreps::Integer = 999,
    seed::Integer = 20260713,
    rng::Union{Nothing, AbstractRNG} = nothing,
)
    nreps >= 2 || throw(ArgumentError("nreps must be at least 2 to calculate sample standard deviations"))
    all(isfinite, tree.edge_length) && all(tree.edge_length .>= 0.0) || throw(ArgumentError("tree contains invalid branch lengths"))
    pool_positions = _nri_tip_positions(tree, collect(pool), "pool")
    observed_positions = _nri_tip_positions(tree, collect(observed), "observed")
    length(observed_positions) >= 2 || throw(ArgumentError("observed must contain at least two tips"))
    length(pool_positions) >= length(observed_positions) || throw(ArgumentError("pool is smaller than observed"))
    pool_set = Set(pool_positions)
    all(in(pool_set), observed_positions) || throw(ArgumentError("observed is not a subset of pool"))

    query_index = _build_tree_query_index(tree)
    observed_mpd, observed_mntd = _community_metrics(tree, observed_positions, query_index)
    mpd_moments = _NriOnlineMoments()
    mntd_moments = _NriOnlineMoments()
    working_pool = copy(pool_positions)
    sampled = Vector{Int}(undef, length(observed_positions))
    swap_positions = Vector{Int}(undef, length(observed_positions))
    random = rng === nothing ? StableRNG(seed) : rng

    for _ in 1:nreps
        for i in eachindex(sampled)
            j = rand(random, i:length(working_pool))
            swap_positions[i] = j
            working_pool[i], working_pool[j] = working_pool[j], working_pool[i]
            sampled[i] = working_pool[i]
        end
        null_mpd, null_mntd = _community_metrics(tree, sampled, query_index)
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
        seed = UInt64(seed),
        status = status,
    )
end
