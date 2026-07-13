Base.@kwdef struct PSVResult
    psv::Vector{Float64} = Float64[]
    richness::Vector{Int} = Int[]
    variance::Vector{Float64} = Float64[]
    scale_vcv::Bool = true
end

Base.@kwdef struct ExpectedPDResult
    richness::Vector{Int} = Int[]
    expected_pd::Vector{Float64} = Float64[]
    variance_pd::Vector{Float64} = Float64[]
    upper_bound::Bool = true
end

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

function _psv_variances(C::AbstractMatrix{<:Real})
    species_count = size(C, 1)
    species_count <= 1 && return fill(NaN, species_count)
    matrix = Matrix{Float64}(C)
    offdiag_mean = sum(matrix) - sum(diag(matrix))
    offdiag_mean /= species_count * (species_count - 1)
    centered = matrix .- offdiag_mean
    for i in 1:species_count
        centered[i, i] = 0.0
    end
    ss1 = sum(abs2, centered) / 2
    ss2 = 0.0
    for i in 1:(species_count - 1), j in (i + 1):species_count
        ss2 += centered[i, j] * (sum(centered[i, :]) - centered[i, j])
    end
    ss3 = -ss1 - ss2
    s1 = 2 * ss1 / (species_count * (species_count - 1))
    s2 = species_count >= 3 ? 2 * ss2 / (species_count * (species_count - 1) * (species_count - 2)) : 0.0
    s3 = species_count >= 4 ? 2 * ss3 / (species_count * (species_count - 1) * (species_count - 2) * (species_count - 3)) : 0.0
    result = fill(NaN, species_count)
    for richness in 2:species_count
        result[richness] = 2 / (richness * (richness - 1)) *
            (s1 + (richness - 2) * s2 + (richness - 2) * (richness - 3) * s3)
    end
    return result
end

function psv(community::AbstractMatrix{<:Real}, tree::CompactTree;
    scale_vcv::Bool = true,
    compute_var::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    abundance = _community_abundance(community, tree, tip_labels)
    C = _phylo_covariance(tree; scale_vcv = scale_vcv)
    variance_by_richness = compute_var ? _psv_variances(C) : fill(NaN, tree.ntips)
    values = fill(NaN, size(abundance, 1))
    richness = zeros(Int, size(abundance, 1))
    variances = fill(NaN, size(abundance, 1))
    for i in axes(abundance, 1)
        present = findall(>(0), abundance[i, :])
        richness[i] = length(present)
        richness[i] > 1 || continue
        sub = C[present, present]
        n = length(present)
        values[i] = (n * sum(diag(sub)) - sum(sub)) / (n * (n - 1))
        compute_var && (variances[i] = variance_by_richness[n])
    end
    return PSVResult(psv = values, richness = richness, variance = variances, scale_vcv = scale_vcv)
end

function psr(community::AbstractMatrix{<:Real}, tree::CompactTree; kwargs...)
    result = psv(community, tree; kwargs...)
    return (psr = result.psv .* result.richness, richness = result.richness, variance = result.variance .* result.richness .^ 2)
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

function edge_abundance_distribution(tree::CompactTree)
    counts = _tools_descendant_tip_counts(tree)
    return (num_children = Int[counts[Int(tree.child_of_edge[e])] for e in 1:tree.nedges], edge_length = copy(tree.edge_length))
end

ead(tree::CompactTree) = edge_abundance_distribution(tree)

function _expected_pd_vectors(tree::CompactTree; upper_bound::Bool = true)
    n = tree.ntips
    desc = _tools_descendant_tip_counts(tree)
    expected = zeros(Float64, n)
    variance = zeros(Float64, n)
    query_index = _build_tree_query_index(tree)
    for richness in 1:n
        p = richness / n
        no_draw = [(1.0 - p)^Int(desc[Int(tree.child_of_edge[edge])]) for edge in 1:tree.nedges]
        for edge in 1:tree.nedges
            k = Int(desc[Int(tree.child_of_edge[edge])])
            q = 1.0 - (1.0 - p)^k
            expected[richness] += tree.edge_length[edge] * q
            variance[richness] += tree.edge_length[edge]^2 * q * (1.0 - q)
        end
        for parent_edge in 1:tree.nedges, descendant_edge in 1:tree.nedges
            parent_edge == descendant_edge && continue
            parent_node = tree.child_of_edge[parent_edge]
            descendant_node = tree.child_of_edge[descendant_edge]
            _tree_is_ancestor(query_index, parent_node, descendant_node) || continue
            covariance = upper_bound ?
                no_draw[parent_edge] * (1.0 - no_draw[parent_edge]) :
                no_draw[parent_edge] * (1.0 - no_draw[descendant_edge])
            variance[richness] += 2.0 * tree.edge_length[parent_edge] * tree.edge_length[descendant_edge] * covariance
        end
    end
    return expected, variance
end

function expected_pd(tree::CompactTree; upper_bound::Bool = true)
    expected, variance = _expected_pd_vectors(tree; upper_bound = upper_bound)
    return ExpectedPDResult(richness = collect(1:tree.ntips), expected_pd = expected, variance_pd = variance, upper_bound = upper_bound)
end

pd_variance(tree::CompactTree; upper_bound::Bool = true) = expected_pd(tree; upper_bound = upper_bound)

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
