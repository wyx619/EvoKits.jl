Base.@kwdef struct PDBootstrapResult
    observed::Vector{Float64} = Float64[]
    replicates::Matrix{Float64} = zeros(Float64, 0, 0)
    richness::Vector{Int} = Int[]
end

function _phylo_covariance(tree::CompactTree; scale_vcv::Bool = true)
    covariance = _tree_shared_path_covariance(tree)
    scale_vcv || return covariance
    diagonal = sqrt.(diag(covariance))
    any(x -> !isfinite(x) || x <= 0, diagonal) && throw(ArgumentError("tree has non-positive root-to-tip variance"))
    return covariance ./ (diagonal * diagonal')
end

function pse(community::AbstractMatrix{<:Real}, tree::CompactTree;
    scale_vcv::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    abundance = _community_abundance(community, tree, tip_labels)
    C = _phylo_covariance(tree; scale_vcv = scale_vcv)
    values = fill(NaN, size(abundance, 1))
    richness = zeros(Int, size(abundance, 1))
    for i in axes(abundance, 1)
        present = findall(>(0), abundance[i, :])
        richness[i] = length(present)
        length(present) > 1 || continue
        weights = abundance[i, present]
        total = sum(weights)
        mean_abundance = mean(weights)
        sub = C[present, present]
        denominator = total^2 - total * mean_abundance
        denominator > 0 || continue
        values[i] = (total * sum(diag(sub) .* weights) - dot(weights, sub * weights)) / denominator
    end
    return (pse = values, richness = richness)
end

function psc(community::AbstractMatrix{<:Real}, tree::CompactTree;
    scale_vcv::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    abundance = _community_abundance(community, tree, tip_labels)
    C = _phylo_covariance(tree; scale_vcv = scale_vcv)
    values = fill(NaN, size(abundance, 1))
    richness = zeros(Int, size(abundance, 1))
    for i in axes(abundance, 1)
        present = findall(>(0), abundance[i, :])
        richness[i] = length(present)
        length(present) > 1 || continue
        sub = Matrix(C[present, present])
        for j in axes(sub, 1)
            sub[j, j] = -1.0
        end
        values[i] = 1.0 - sum(maximum(sub, dims = 2)) / length(present)
    end
    return (psc = values, richness = richness)
end

function psd(community::AbstractMatrix{<:Real}, tree;
    compute_var::Bool = true,
    scale_vcv::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    v = psv(community, tree; compute_var = compute_var, scale_vcv = scale_vcv, tip_labels = tip_labels)
    r = psr(community, tree; compute_var = compute_var, scale_vcv = scale_vcv, tip_labels = tip_labels)
    e = pse(community, tree; scale_vcv = scale_vcv, tip_labels = tip_labels)
    c = psc(community, tree; scale_vcv = scale_vcv, tip_labels = tip_labels)
    return (psv = v.psv, var_psv = v.variance, psc = c.psc, psr = r.psr, var_psr = r.variance, pse = e.pse, richness = v.richness)
end

function _normalize_pd_variant(variant)
    value = Symbol(lowercase(String(variant)))
    value in (:tbl, :mst, :ueh, :sbl, :tip) || throw(ArgumentError("variant must be TBL, MST, UEH, SBL, or TIP"))
    return value
end

function pd_calc(community::AbstractMatrix{<:Real}, tree::CompactTree;
    variant = :tbl,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    variant = _normalize_pd_variant(variant)
    abundance = _community_abundance(community, tree, tip_labels)
    active = _community_edge_presence(abundance, tree)
    values = zeros(Float64, size(abundance, 1))
    terminal = Set(Int32.(tree.tip_ids))
    query_index = _build_tree_query_index(tree)
    for i in axes(abundance, 1)
        present = findall(>(0), abundance[i, :])
        present_nodes = Int32[tree.tip_ids[j] for j in present]
        mrca_node = isempty(present_nodes) ? Int32(0) : foldl((a, b) -> _tree_lca(query_index, a, b), present_nodes)
        counts = zeros(Int, tree.nedges)
        for edge in 1:tree.nedges
            child = tree.child_of_edge[edge]
            counts[edge] = count(tip -> _tree_is_ancestor(query_index, child, tree.tip_ids[tip]), present)
        end
        for edge in 1:tree.nedges
            child = tree.child_of_edge[edge]
            include = if variant === :tbl
                active[i, edge]
            elseif variant === :sbl
                active[i, edge] && counts[edge] > 1
            elseif variant === :ueh
                active[i, edge] && counts[edge] == 1
            elseif variant === :tip
                active[i, edge] && child in terminal
            else
                active[i, edge] && mrca_node != 0 && child != mrca_node && _tree_is_ancestor(query_index, mrca_node, child)
            end
            if include
                values[i] += tree.edge_length[edge]
            end
        end
    end
    return (pd = values, richness = vec(sum(abundance .> 0, dims = 2)), variant = variant)
end

function pd_calc(tree::CompactTree, tips;
    variant = :tbl,
)
    row = zeros(Float64, tree.ntips)
    values = tips isa AbstractString ? [String(tips)] : collect(tips)
    length(unique(values)) == length(values) || throw(ArgumentError("tips must not contain duplicates"))
    for value in values
        position = if value isa Integer
            1 <= value <= tree.ntips || throw(ArgumentError("tip position is outside the tree"))
            Int(value)
        else
            node = get(tree.tipname_to_id, String(value), Int32(0))
            node == 0 && throw(ArgumentError("unknown tree tip: $value"))
            findfirst(==(node), tree.tip_ids)
        end
        row[position] = 1.0
    end
    result = pd_calc(reshape(row, 1, :), tree; variant = variant)
    return result.pd[1]
end

function pd_bootstrap(community::AbstractMatrix{<:Real}, tree::CompactTree;
    reps::Integer = 999,
    tip_weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
    include_root::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    rng::AbstractRNG = Random.default_rng(),
)
    reps >= 1 || throw(ArgumentError("reps must be positive"))
    abundance = _community_abundance(community, tree, tip_labels)
    weights = tip_weights === nothing ? ones(Float64, tree.ntips) : Float64.(tip_weights)
    length(weights) == tree.ntips || throw(ArgumentError("tip_weights must match tree tips"))
    all(isfinite, weights) && all(>=(0), weights) && sum(weights) > 0 || throw(ArgumentError("tip_weights must be finite and non-negative with positive sum"))
    probabilities = weights / sum(weights)
    observed = pd(abundance, tree; include_root = include_root).pd
    draws = Matrix{Float64}(undef, size(abundance, 1), reps)
    for r in 1:reps
        sampled = falses(size(abundance))
        for i in axes(abundance, 1)
            richness = count(>(0), abundance[i, :])
            richness == 0 && continue
            keys = [probabilities[j] > 0 ? log(rand(rng)) / probabilities[j] : -Inf for j in 1:tree.ntips]
            order = sortperm(keys; rev = true)
            sampled[i, order[1:richness]] .= true
        end
        draws[:, r] = pd(Float64.(sampled), tree; include_root = include_root).pd
    end
    return PDBootstrapResult(observed = observed, replicates = draws, richness = vec(sum(abundance .> 0, dims = 2)))
end

function pd_bootstrap(tree::CompactTree, ntips::Integer;
    reps::Integer = 1000,
    method = :TBL,
    tip_weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
    rng::AbstractRNG = Random.default_rng(),
)
    1 <= ntips <= tree.ntips || throw(ArgumentError("ntips must be between 1 and tree.ntips"))
    reps >= 1 || throw(ArgumentError("reps must be positive"))
    weights = tip_weights === nothing ? ones(Float64, tree.ntips) : Float64.(tip_weights)
    length(weights) == tree.ntips || throw(ArgumentError("tip_weights must match tree tips"))
    all(isfinite, weights) && all(>=(0), weights) && sum(weights) > 0 || throw(ArgumentError("tip_weights must be finite and non-negative with positive sum"))
    probabilities = weights / sum(weights)
    values = Vector{Float64}(undef, reps)
    for r in 1:reps
        keys = [probabilities[j] > 0 ? log(rand(rng)) / probabilities[j] : -Inf for j in 1:tree.ntips]
        tips = sortperm(keys; rev = true)[1:ntips]
        values[r] = pd_calc(tree, tips; variant = method)
    end
    return values
end
