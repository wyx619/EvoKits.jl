Base.@kwdef struct PICResult
    contrasts::Vector{Float64} = Float64[]
    nodes::Vector{Int32} = Int32[]
    distances::Vector{Float64} = Float64[]
    ancestral_values::Vector{Float64} = Float64[]
    tip_labels::Vector{String} = String[]
    scaled::Bool = true
    weighted::Bool = true
end

function phylogenetic_contrasts(tree::CompactTree, trait::AbstractVector{<:Real};
    scaled::Bool = true,
    weighted::Bool = true,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    values = _trait_in_tree_order(tree, trait, tip_labels)
    tree.ntips >= 2 || throw(ArgumentError("PIC requires at least two tips"))
    all(isfinite, tree.edge_length) && all(>(0), tree.edge_length) ||
        throw(ArgumentError("PIC requires finite, strictly positive branch lengths"))
    for node in tree.postorder_internal
        length(tree.children[Int(node)]) == 2 || throw(ArgumentError("PIC currently requires a strictly bifurcating tree"))
    end

    edge_of_child = _tip_edge_index(tree)
    state = fill(NaN, tree.nnodes)
    variance = zeros(Float64, tree.nnodes)
    for (i, tip) in enumerate(tree.tip_ids)
        state[Int(tip)] = values[i]
        variance[Int(tip)] = 0.0
    end
    contrasts = Float64[]
    nodes = Int32[]
    distances = Float64[]
    ancestral = Float64[]
    for node32 in tree.postorder_internal
        node = Int(node32)
        child1, child2 = Int.(tree.children[node])
        edge1 = Int(edge_of_child[child1])
        edge2 = Int(edge_of_child[child2])
        length1 = weighted ? tree.edge_length[edge1] : 1.0
        length2 = weighted ? tree.edge_length[edge2] : 1.0
        v1 = variance[child1] + length1
        v2 = variance[child2] + length2
        total = v1 + v2
        total > 0 || throw(ArgumentError("PIC encountered a zero contrast variance"))
        raw = state[child1] - state[child2]
        contrast = scaled ? raw / sqrt(total) : raw
        push!(contrasts, contrast)
        push!(nodes, node32)
        push!(distances, total)
        state[node] = (state[child1] * v2 + state[child2] * v1) / total
        variance[node] = (v1 * v2) / total
        push!(ancestral, state[node])
    end
    return PICResult(
        contrasts = contrasts,
        nodes = nodes,
        distances = distances,
        ancestral_values = ancestral,
        tip_labels = copy(tree.tip_labels),
        scaled = scaled,
        weighted = weighted,
    )
end
