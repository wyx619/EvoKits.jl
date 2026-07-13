Base.@kwdef struct PhyloSignalResult
    k::Float64 = NaN
    pic_variance_observed::Float64 = NaN
    pic_variance_null_mean::Float64 = NaN
    pic_variance_null_sd::Float64 = NaN
    p_value::Float64 = NaN
    z_score::Float64 = NaN
    trait::Vector{Float64} = Float64[]
    tip_labels::Vector{String} = String[]
    valid_replicates::Int = 0
    permutations::Int = 0
    status::Symbol = :ok
end

function _trait_in_tree_order(tree::CompactTree, trait::AbstractVector{<:Real}, tip_labels)
    length(trait) == tree.ntips || throw(ArgumentError("trait length must match tree tips"))
    values = Float64.(trait)
    all(isfinite, values) || throw(ArgumentError("trait contains non-finite values"))
    tip_labels === nothing && return values
    length(tip_labels) == tree.ntips || throw(ArgumentError("tip_labels length must match trait length"))
    positions = Dict{String, Int}(String(label) => i for (i, label) in enumerate(tip_labels))
    all(haskey(positions, label) for label in tree.tip_labels) || throw(ArgumentError("tip_labels must contain every tree tip exactly once"))
    length(positions) == tree.ntips || throw(ArgumentError("tip_labels must not contain duplicates"))
    return [values[positions[label]] for label in tree.tip_labels]
end

function _pic_variance(tree::CompactTree, trait::AbstractVector{<:Real})
    result = phylogenetic_contrasts(tree, trait; scaled = true, weighted = true)
    n = length(result.contrasts)
    n > 1 || return NaN
    return sum(abs2, result.contrasts) / (n - 1)
end

function _phylo_k(tree::CompactTree, trait::Vector{Float64})
    n = tree.ntips
    n > 2 || return NaN
    C = _tree_shared_path_covariance(tree)
    precision = pinv(Symmetric(C))
    denominator = sum(precision)
    denominator > 0 || return NaN
    mean_phylo = sum(precision * trait) / denominator
    residual = trait .- mean_phylo
    mse = dot(residual, precision * residual) / (n - 1)
    mse0 = sum((trait .- mean(trait)) .^ 2) / (n - 1)
    expected = (sum(diag(C)) - n / denominator) / (n - 1)
    (mse > 0 && expected > 0) || return NaN
    return (mse0 / mse) / expected
end

function phylogenetic_signal_k(tree::CompactTree, trait::AbstractVector{<:Real};
    permutations::Integer = 999,
    rng::AbstractRNG = StableRNG(20260713),
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    permutations >= 0 || throw(ArgumentError("permutations must be non-negative"))
    values = _trait_in_tree_order(tree, trait, tip_labels)
    observed_k = _phylo_k(tree, values)
    observed_variance = _pic_variance(tree, values)
    null = Float64[]
    working = similar(values)
    for _ in 1:permutations
        copyto!(working, values)
        Random.shuffle!(rng, working)
        value = _pic_variance(tree, working)
        isfinite(value) && push!(null, value)
    end
    null_mean = isempty(null) ? NaN : mean(null)
    null_sd = length(null) > 1 ? std(null) : NaN
    rank = if isfinite(observed_variance) && !isempty(null)
        1.0 + count(x -> x < observed_variance, null) + 0.5 * count(==(observed_variance), null)
    else
        NaN
    end
    p = isfinite(rank) ? rank / (length(null) + 1) : NaN
    z = isfinite(null_sd) && null_sd > 0 ? (observed_variance - null_mean) / null_sd : NaN
    return PhyloSignalResult(
        k = observed_k,
        pic_variance_observed = observed_variance,
        pic_variance_null_mean = null_mean,
        pic_variance_null_sd = null_sd,
        p_value = p,
        z_score = z,
        trait = values,
        tip_labels = copy(tree.tip_labels),
        valid_replicates = length(null),
        permutations = Int(permutations),
        status = isempty(null) && permutations > 0 ? :no_valid_replicates : :ok,
    )
end
