@inline function _tools_sample_labels(n::Integer, labels::Union{Nothing, AbstractVector{<:AbstractString}})
    labels === nothing && return ["sample_$i" for i in 1:n]
    length(labels) == n || throw(ArgumentError("sample_labels length must match number of sample rows"))
    return String.(labels)
end

function _tools_tip_permutation(tree::CompactTree, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}})
    tip_labels === nothing && return collect(1:tree.ntips)
    length(tip_labels) == tree.ntips || throw(ArgumentError("tip_labels length must match number of sample columns"))
    labels = String.(tip_labels)
    length(unique(labels)) == length(labels) || throw(ArgumentError("tip_labels must not contain duplicates"))
    perm = Vector{Int}(undef, tree.ntips)
    for (j, label) in enumerate(labels)
        node = get(tree.tipname_to_id, label, Int32(0))
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
    all(x -> isfinite(x) && x >= 0, samp) || throw(ArgumentError("community values must be finite and non-negative"))
    abund = Matrix{Float64}(undef, size(samp, 1), length(perm))
    @inbounds for i in axes(samp, 1), j in axes(samp, 2)
        abund[i, perm[j]] = Float64(samp[i, j])
    end
    return abund
end

function _tools_presence_only(samp::AbstractMatrix{<:Real}, perm::Vector{Int})
    all(x -> isfinite(x) && x >= 0, samp) || throw(ArgumentError("community values must be finite and non-negative"))
    presence = Matrix{Bool}(undef, size(samp, 1), length(perm))
    @inbounds for i in axes(samp, 1), j in axes(samp, 2)
        presence[i, perm[j]] = samp[i, j] > 0
    end
    return presence
end

function _tools_presence_matrix(samp::AbstractMatrix{<:Real}, perm::Vector{Int})
    return _tools_presence_only(samp, perm)
end

function _community_abundance(samp::AbstractMatrix{<:Real}, tree::CompactTree, tip_labels)
    size(samp, 2) == tree.ntips || throw(ArgumentError("community matrix must have one column per tree tip"))
    perm = _tools_tip_permutation(tree, tip_labels)
    return _tools_abundance_only(samp, perm)
end

function _community_presence(samp::AbstractMatrix{<:Real}, tree::CompactTree, tip_labels)
    abundance = _community_abundance(samp, tree, tip_labels)
    return abundance .> 0
end
