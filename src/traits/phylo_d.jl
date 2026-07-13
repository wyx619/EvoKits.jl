"""
    PhyloDResult

Result of a Fritz-Purvis phylogenetic D calculation for a binary trait.
The random and Brownian vectors contain one sum of nodal changes per
permutation. Nodal matrices are populated when `store_nodal_values=true`.
"""
Base.@kwdef struct PhyloDResult
    d::Float64 = NaN
    p_random::Float64 = NaN
    p_brownian::Float64 = NaN
    observed::Float64 = NaN
    mean_random::Float64 = NaN
    mean_brownian::Float64 = NaN
    state_levels::Vector{String} = String[]
    state_counts::Vector{Int} = Int[]
    tip_labels::Vector{String} = String[]
    node_ids::Vector{Int32} = Int32[]
    observed_nodal_values::Vector{Float64} = Float64[]
    random_changes::Vector{Float64} = Float64[]
    brownian_changes::Vector{Float64} = Float64[]
    random_nodal_values::Union{Nothing, Matrix{Float64}} = nothing
    brownian_nodal_values::Union{Nothing, Matrix{Float64}} = nothing
    permutations::Int = 0
    random_weights::Union{Nothing, Vector{Float64}} = nothing
    status::Symbol = :ok
end

# Keep the names used by caper available for callers migrating old workflows.
function Base.getproperty(result::PhyloDResult, name::Symbol)
    name === :D && return getfield(result, :d)
    name === :DEstimate && return getfield(result, :d)
    name === :Pval1 && return getfield(result, :p_random)
    name === :Pval0 && return getfield(result, :p_brownian)
    name === :Parameters && return (
        Observed = getfield(result, :observed),
        MeanRandom = getfield(result, :mean_random),
        MeanBrownian = getfield(result, :mean_brownian),
    )
    name === :Permutations && return (
        random = getfield(result, :random_changes),
        brownian = getfield(result, :brownian_changes),
    )
    name === :NodalVals && return (
        observed = getfield(result, :observed_nodal_values),
        random = getfield(result, :random_nodal_values),
        brownian = getfield(result, :brownian_nodal_values),
    )
    return getfield(result, name)
end

function _phylo_d_binary_trait(trait)
    values = collect(trait)
    isempty(values) && throw(ArgumentError("binary trait cannot be empty"))
    any(ismissing, values) && throw(ArgumentError("binary trait contains missing values"))

    if all(value -> value isa Real, values)
        numeric = Float64.(values)
        all(isfinite, numeric) || throw(ArgumentError("binary trait contains non-finite values"))
        levels = sort(unique(numeric))
        length(levels) == 2 || throw(ArgumentError("binary trait must contain exactly two states"))
        state = Int.(numeric .== levels[2])
        return state, string.(levels)
    end

    labels = String.(values)
    levels = sort(unique(labels))
    length(levels) == 2 || throw(ArgumentError("binary trait must contain exactly two states"))
    state = Int.(labels .== levels[2])
    return state, levels
end

function _phylo_d_validate_tree(tree::CompactTree)
    tree.ntips >= 2 || throw(ArgumentError("phylo_d requires at least two tree tips"))
    all(isfinite, tree.edge_length) && all(tree.edge_length .> 0.0) ||
        throw(ArgumentError("phylo_d requires finite, strictly positive branch lengths"))
    for node in tree.postorder_internal
        length(tree.children[Int(node)]) >= 2 ||
            throw(ArgumentError("phylo_d does not support unary internal nodes"))
    end
    return nothing
end

function _phylo_d_change!(tree::CompactTree, tip_states::AbstractVector{<:Real}, edge_to_node::Vector{Float64}, node_values::Vector{Float64}, nodal_values::Vector{Float64}, internal_index::Vector{Int})
    length(tip_states) == tree.ntips || throw(ArgumentError("tip state length does not match tree"))
    @inbounds for (tip_position, tip) in enumerate(tree.tip_ids)
        node_values[Int(tip)] = Float64(tip_states[tip_position])
    end

    total = 0.0
    @inbounds for node32 in tree.postorder_internal
        node = Int(node32)
        children = tree.children[node]
        nchildren = length(children)
        branch_lengths = Vector{Float64}(undef, nchildren)
        child_values = Vector{Float64}(undef, nchildren)
        for i in 1:nchildren
            child = Int(children[i])
            branch_lengths[i] = edge_to_node[child]
            child_values[i] = node_values[child]
        end

        invbl = 1.0 ./ branch_lengths
        if nchildren == 2
            weight_sum = sum(invbl)
            node_value = sum(child_values .* invbl) / weight_sum
            branch_adjustment = prod(branch_lengths) / sum(branch_lengths)
        else
            # This is the polytomy treatment used by caper's crunch engine.
            group = if var(child_values) < eps(Float64)
                [true; falses(nchildren - 1)]
            else
                child_values .> mean(child_values)
            end
            any(group) && any(.!group) || throw(ArgumentError("could not split a polytomy into two state groups"))

            sub_weights = [sum(invbl[.!group]), sum(invbl[group])]
            sub_values = [
                sum(child_values[.!group] .* invbl[.!group]) / sub_weights[1],
                sum(child_values[group] .* invbl[group]) / sub_weights[2],
            ]
            node_value = sum(sub_values .* sub_weights) / sum(sub_weights)
            branch_adjustment = 1.0 / sum(sub_weights)
        end

        node_values[node] = node_value
        nodal_values[internal_index[node]] = node_value
        total += sum(abs.(child_values .- node_value))
        node != Int(tree.root) && (edge_to_node[node] += branch_adjustment)
    end
    return total
end

function _phylo_d_covariance(tree::CompactTree)
    index = _build_tree_query_index(tree)
    covariance = Matrix{Float64}(undef, tree.ntips, tree.ntips)
    @inbounds for i in 1:tree.ntips
        first_tip = tree.tip_ids[i]
        for j in 1:i
            second_tip = tree.tip_ids[j]
            shared_node = _tree_lca(index, first_tip, second_tip)
            value = tree.dist_from_root[Int(shared_node)]
            covariance[i, j] = value
            covariance[j, i] = value
        end
    end
    return covariance
end

function _phylo_d_random_states(rng::AbstractRNG, states::Vector{Int}, weights::Union{Nothing, Vector{Float64}})
    permutation = if weights === nothing
        randperm(rng, length(states))
    else
        keys = Vector{Float64}(undef, length(states))
        @inbounds for i in eachindex(weights)
            keys[i] = weights[i] > 0.0 ? log(rand(rng)) / weights[i] : -Inf
        end
        sortperm(keys; rev = true)
    end
    return states[permutation]
end

"""
    phylo_d(tree, trait; permutations=1000, rng=StableRNG(20260713),
            random_weights=nothing, store_nodal_values=false)

Compute the Fritz-Purvis D statistic for a binary trait. The random null model
shuffles observed tip states without replacement. The Brownian threshold null
simulates a multivariate normal trait using shared root-to-tip branch lengths,
then thresholds each replicate to the observed state prevalence.
"""
function phylo_d(
    tree::CompactTree,
    trait;
    permutations::Integer = 1000,
    rng::AbstractRNG = StableRNG(20260713),
    random_weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
    store_nodal_values::Bool = false,
)
    permutations >= 1 || throw(ArgumentError("permutations must be positive"))
    _phylo_d_validate_tree(tree)
    state, state_levels = _phylo_d_binary_trait(trait)
    length(state) == tree.ntips || throw(ArgumentError("binary trait length must match tree.ntips"))

    weights = if random_weights === nothing
        nothing
    else
        values = Float64.(random_weights)
        length(values) == tree.ntips || throw(ArgumentError("random_weights length must match tree.ntips"))
        all(isfinite, values) && all(values .>= 0.0) && sum(values) > 0.0 ||
            throw(ArgumentError("random_weights must be finite, non-negative, and not all zero"))
        values
    end

    internal_nodes = copy(tree.postorder_internal)
    internal_index = zeros(Int, tree.nnodes)
    @inbounds for (i, node) in enumerate(internal_nodes)
        internal_index[Int(node)] = i
    end
    ninternal = length(internal_nodes)
    edge_lengths = zeros(Float64, tree.nnodes)
    @inbounds for edge in 1:tree.nedges
        edge_lengths[Int(tree.child_of_edge[edge])] = tree.edge_length[edge]
    end
    node_values = zeros(Float64, tree.nnodes)
    nodal_work = zeros(Float64, ninternal)

    observed_edges = copy(edge_lengths)
    observed_nodal = zeros(Float64, ninternal)
    observed_changes = _phylo_d_change!(tree, state, observed_edges, node_values, observed_nodal, internal_index)

    random_changes = Vector{Float64}(undef, Int(permutations))
    brownian_changes = Vector{Float64}(undef, Int(permutations))
    random_nodal = store_nodal_values ? Matrix{Float64}(undef, ninternal, Int(permutations)) : nothing
    brownian_nodal = store_nodal_values ? Matrix{Float64}(undef, ninternal, Int(permutations)) : nothing

    state_prevalence = count(==(1), state) / length(state)
    for replicate in 1:Int(permutations)
        random_state = _phylo_d_random_states(rng, state, weights)
        random_edges = copy(edge_lengths)
        random_changes[replicate] = _phylo_d_change!(tree, random_state, random_edges, node_values, nodal_work, internal_index)
        store_nodal_values && (random_nodal[:, replicate] .= nodal_work)
    end

    covariance = _phylo_d_covariance(tree)
    brownian_distribution = MvNormal(zeros(Float64, tree.ntips), Symmetric(covariance))
    liabilities = rand(rng, brownian_distribution, Int(permutations))
    for replicate in 1:Int(permutations)
        threshold = quantile(@view(liabilities[:, replicate]), state_prevalence)
        brownian_state = Int.(@view(liabilities[:, replicate]) .< threshold)
        brownian_edges = copy(edge_lengths)
        brownian_changes[replicate] = _phylo_d_change!(tree, brownian_state, brownian_edges, node_values, nodal_work, internal_index)
        store_nodal_values && (brownian_nodal[:, replicate] .= nodal_work)
    end

    mean_random = mean(random_changes)
    mean_brownian = mean(brownian_changes)
    denominator = mean_random - mean_brownian
    d = abs(denominator) > eps(Float64) ? (observed_changes - mean_brownian) / denominator : NaN
    p_random = count(value -> value < observed_changes, random_changes) / Float64(permutations)
    p_brownian = count(value -> value > observed_changes, brownian_changes) / Float64(permutations)
    status = isfinite(d) ? :ok : :null_range_zero

    return PhyloDResult(
        d = d,
        p_random = p_random,
        p_brownian = p_brownian,
        observed = observed_changes,
        mean_random = mean_random,
        mean_brownian = mean_brownian,
        state_levels = state_levels,
        state_counts = [count(==(0), state), count(==(1), state)],
        tip_labels = copy(tree.tip_labels),
        node_ids = internal_nodes,
        observed_nodal_values = observed_nodal,
        random_changes = random_changes,
        brownian_changes = brownian_changes,
        random_nodal_values = random_nodal,
        brownian_nodal_values = brownian_nodal,
        permutations = Int(permutations),
        random_weights = weights,
        status = status,
    )
end

function phylo_d(
    tree::CompactTree,
    data::AbstractDataFrame,
    trait_col;
    taxon_col::Union{Nothing, Symbol, AbstractString} = nothing,
    kwargs...,
)
    trait = align_traits_to_tree(tree, data; taxon_col = taxon_col, trait_cols = trait_col)
    return phylo_d(tree, trait; kwargs...)
end

"""
    PhyloDSubsetResult

Phylogenetic D results calculated independently for filtered internal clades.
Entries in `raw` are `nothing` when a clade contains only one binary state.
"""
Base.@kwdef struct PhyloDSubsetResult
    raw::Vector{Union{Nothing, PhyloDResult}} = Union{Nothing, PhyloDResult}[]
    d::Vector{Float64} = Float64[]
    p_random::Vector{Float64} = Float64[]
    p_brownian::Vector{Float64} = Float64[]
    phylogenetic_depth::Vector{Float64} = Float64[]
    tips::Vector{Int} = Int[]
    nodes::Vector{Int} = Int[]
    state_count::Vector{Int} = Int[]
    tip_sets::Vector{Vector{String}} = Vector{Vector{String}}()
end

function Base.getproperty(result::PhyloDSubsetResult, name::Symbol)
    name === :D && return getfield(result, :d)
    name === :DEstimate && return getfield(result, :d)
    name === :Pval1 && return getfield(result, :p_random)
    name === :Pval0 && return getfield(result, :p_brownian)
    return getfield(result, name)
end

"""
    phylo_d_subset(tree, trait; permutations=1000, min_tips=1,
                   max_tips=tree.ntips, min_nodes=1,
                   max_nodes=length(tree.postorder_internal), kwargs...)

Calculate `phylo_d` for each internal clade passing the tip and node filters.
Clades with only one observed state are retained in the output with `NaN`
summary values and `nothing` in `raw`, matching `caper::phylo.d.subset`.
"""
function phylo_d_subset(
    tree::CompactTree,
    trait;
    permutations::Integer = 1000,
    min_tips::Integer = 1,
    max_tips::Integer = tree.ntips,
    min_nodes::Integer = 1,
    max_nodes::Integer = length(tree.postorder_internal),
    rng::AbstractRNG = StableRNG(20260713),
    random_weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
    store_nodal_values::Bool = false,
)
    min_tips >= 1 && max_tips >= min_tips || throw(ArgumentError("invalid tip-count bounds"))
    min_nodes >= 1 && max_nodes >= min_nodes || throw(ArgumentError("invalid node-count bounds"))
    _phylo_d_validate_tree(tree)
    state, _ = _phylo_d_binary_trait(trait)
    length(state) == tree.ntips || throw(ArgumentError("binary trait length must match tree.ntips"))

    tip_position_of_node = zeros(Int, tree.nnodes)
    @inbounds for (position, tip) in enumerate(tree.tip_ids)
        tip_position_of_node[Int(tip)] = position
    end
    descendant_tips = [Int[] for _ in 1:tree.nnodes]
    descendant_nodes = zeros(Int, tree.nnodes)
    @inbounds for node32 in tree.postorder
        node = Int(node32)
        if tree.is_tip[node]
            push!(descendant_tips[node], tip_position_of_node[node])
        else
            for child32 in tree.children[node]
                child = Int(child32)
                append!(descendant_tips[node], descendant_tips[child])
                descendant_nodes[node] += descendant_nodes[child]
            end
            descendant_nodes[node] += 1
        end
    end

    candidates = Int32[]
    for node32 in tree.postorder_internal
        node = Int(node32)
        ntips = length(descendant_tips[node])
        nnodes = descendant_nodes[node]
        min_tips <= ntips <= max_tips && min_nodes <= nnodes <= max_nodes && push!(candidates, node32)
    end

    raw = Vector{Union{Nothing, PhyloDResult}}(undef, length(candidates))
    d = fill(NaN, length(candidates))
    p_random = fill(NaN, length(candidates))
    p_brownian = fill(NaN, length(candidates))
    phylogenetic_depth = fill(NaN, length(candidates))
    tips = zeros(Int, length(candidates))
    nodes = zeros(Int, length(candidates))
    state_count = zeros(Int, length(candidates))
    tip_sets = Vector{Vector{String}}(undef, length(candidates))

    for (i, node32) in enumerate(candidates)
        node = Int(node32)
        positions = descendant_tips[node]
        labels = tree.tip_labels[positions]
        subtree = keep_tip(tree, labels)
        subtree_positions = [findfirst(==(label), tree.tip_labels) for label in subtree.tip_labels]
        substate = state[subtree_positions]
        tip_sets[i] = copy(labels)
        tips[i] = length(positions)
        nodes[i] = descendant_nodes[node]
        state_count[i] = count(==(1), substate)
        phylogenetic_depth[i] = subtree.dist_from_root[Int(subtree.tip_ids[1])]

        if length(unique(substate)) < 2
            raw[i] = nothing
            continue
        end

        subweights = if random_weights === nothing
            nothing
        else
            all_weights = Float64.(random_weights)
            length(all_weights) == tree.ntips || throw(ArgumentError("random_weights length must match tree.ntips"))
            all_weights[subtree_positions]
        end
        raw[i] = phylo_d(
            subtree,
            substate;
            permutations = permutations,
            rng = rng,
            random_weights = subweights,
            store_nodal_values = store_nodal_values,
        )
        d[i] = raw[i].d
        p_random[i] = raw[i].p_random
        p_brownian[i] = raw[i].p_brownian
    end

    return PhyloDSubsetResult(
        raw = raw,
        d = d,
        p_random = p_random,
        p_brownian = p_brownian,
        phylogenetic_depth = phylogenetic_depth,
        tips = tips,
        nodes = nodes,
        state_count = state_count,
        tip_sets = tip_sets,
    )
end

function phylo_d_subset(
    tree::CompactTree,
    data::AbstractDataFrame,
    trait_col;
    taxon_col::Union{Nothing, Symbol, AbstractString} = nothing,
    kwargs...,
)
    trait = align_traits_to_tree(tree, data; taxon_col = taxon_col, trait_cols = trait_col)
    return phylo_d_subset(tree, trait; kwargs...)
end
