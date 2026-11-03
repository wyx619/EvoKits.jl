Base.@kwdef struct EvolutionaryDistinctivenessResult
    tip_labels::Vector{String} = String[]
    scores::Vector{Float64} = Float64[]
    type::Symbol = :unknown
    scale::Bool = false
    use_branch_lengths::Bool = true
    polytomy_correction::Symbol = :none
end

@inline function _ed_normalize_type(type::Symbol)
    type === :equal_splits && return type
    type === :fair_proportion && return type
    type === :fair_proportions && return :fair_proportion
    type === :fairproportion && return :fair_proportion
    type === :fp && return :fair_proportion
    type === :es && return :equal_splits
    throw(ArgumentError("type must be :equal_splits or :fair_proportion"))
end

@inline function _ed_scale_factor(tree::CompactTree)
    total = 0.0
    @inbounds for edge in 1:tree.nedges
        len = tree.edge_length[edge]
        isfinite(len) || throw(ArgumentError("tree contains non-finite branch lengths"))
        total += len
    end
    total > 0.0 || throw(ArgumentError("cannot scale a tree with zero total branch length"))
    tip_depth = 0.0
    ultrametric = true
    if tree.ntips > 0
        first_tip = Int(tree.tip_ids[1])
        tip_depth = tree.dist_from_root[first_tip]
        atol = max(1e-8, 1e-8 * max(abs(tip_depth), 1.0))
        @inbounds for tip in tree.tip_ids
            if abs(tree.dist_from_root[Int(tip)] - tip_depth) > atol
                ultrametric = false
                break
            end
        end
    end
    if ultrametric
        tip_depth > 0.0 || throw(ArgumentError("cannot scale an ultrametric tree with zero depth"))
        return 1.0 / tip_depth
    end
    return 1.0 / total
end

function _ed_normalize_polytomy_correction(value)
    correction = Symbol(lowercase(String(value)))
    correction in (:isaac, :mooers, :none) || throw(ArgumentError("polytomy_correction must be :isaac, :mooers, or :none"))
    return correction
end

function _ed_correction_factor(parent_degree::Int, correction::Symbol)
    parent_degree <= 2 && return 1.0
    correction === :none && return 1.0
    if correction === :isaac
        parent_degree > 57 && return 0.0
        return 1.081 - 0.267 * log(parent_degree)
    end
    harmonic = sum(inv(k) for k in 2:parent_degree)
    return (parent_degree - 1) / (parent_degree * harmonic)
end

function _ed_score_fair_proportion!(scores::Vector{Float64}, tree::CompactTree, desc_count::Vector{Int32}, scale_factor::Float64, use_branch_lengths::Bool, polytomy_correction::Symbol)
    acc = zeros(Float64, tree.nnodes)
    @inbounds for node in tree.preorder
        inode = Int(node)
        base = acc[inode]
        first_edge = Int(tree.first_child_edge[inode])
        first_edge == 0 && continue
        last_edge = Int(tree.last_child_edge[inode])
        for edge in first_edge:last_edge
            child = Int(tree.child_of_edge[edge])
            len = use_branch_lengths ? tree.edge_length[edge] * scale_factor : 1.0
            len *= _ed_correction_factor(last_edge - first_edge + 1, polytomy_correction)
            acc[child] = base + len / Float64(desc_count[child])
        end
    end
    @inbounds for (i, tip) in enumerate(tree.tip_ids)
        scores[i] = acc[Int(tip)]
    end
    return scores
end

function _ed_score_equal_splits!(scores::Vector{Float64}, tree::CompactTree, scale_factor::Float64, use_branch_lengths::Bool)
    # Picante assigns the smallest 0.5^n weight to the root-most internal
    # edge and the largest 0.5 weight to the internal edge nearest the tip.
    # Walk each tip path upward so this ordering is explicit and independent
    # of the traversal used to build CompactTree.
    edge_of_child = zeros(Int, tree.nnodes)
    @inbounds for edge in 1:tree.nedges
        edge_of_child[Int(tree.child_of_edge[edge])] = edge
    end

    @inbounds for (i, tip_id) in enumerate(tree.tip_ids)
        tip = Int(tip_id)
        pendant_edge = edge_of_child[tip]
        score = use_branch_lengths ? tree.edge_length[pendant_edge] * scale_factor : 1.0
        node = Int(tree.parent_of_node[tip])
        weight = 0.5
        while node != 0 && node != Int(tree.root)
            edge = edge_of_child[node]
            len = use_branch_lengths ? tree.edge_length[edge] * scale_factor : 1.0
            score += len * weight
            weight *= 0.5
            node = Int(tree.parent_of_node[node])
        end
        scores[i] = score
    end
    return scores
end

function evol_distinct(tree::CompactTree; type::Symbol = :equal_splits, scale::Bool = false, use_branch_lengths::Bool = true, polytomy_correction = :none)
    ed_type = _ed_normalize_type(type)
    correction = _ed_normalize_polytomy_correction(polytomy_correction)
    tree.ntips > 0 || throw(ArgumentError("tree must contain at least one tip"))
    scale_factor = (scale && use_branch_lengths) ? _ed_scale_factor(tree) : 1.0
    scores = Vector{Float64}(undef, tree.ntips)
    if ed_type === :fair_proportion
        desc_count = _tools_descendant_tip_counts(tree)
        _ed_score_fair_proportion!(scores, tree, desc_count, scale_factor, use_branch_lengths, correction)
    else
        _ed_score_equal_splits!(scores, tree, scale_factor, use_branch_lengths)
    end
    return EvolutionaryDistinctivenessResult(
        tip_labels = copy(tree.tip_labels),
        scores = scores,
        type = ed_type,
        scale = scale,
        use_branch_lengths = use_branch_lengths,
        polytomy_correction = correction,
    )
end

function evol_distinct_dataframe(tree::CompactTree; kwargs...)
    result = evol_distinct(tree; kwargs...)
    return DataFrame(Species = result.tip_labels, w = result.scores)
end
