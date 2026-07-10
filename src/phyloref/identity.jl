function _phyloref_first_tip_label(tree::CompactTree, map::PhyloRef, node::Integer)
    1 <= node <= tree.nnodes || throw(ArgumentError("node $node is outside 1:$(tree.nnodes)"))
    positions = map.descendant_tip_positions_by_node[Int(node)]
    isempty(positions) && throw(ArgumentError("node $node has no descendant tips"))
    return tree.tip_labels[Int(first(positions))]
end

function phylo_node_anchor(tree::CompactTree, node::Integer; map::PhyloRef = build_phyloref(tree))
    1 <= node <= tree.nnodes || throw(ArgumentError("node $node is outside 1:$(tree.nnodes)"))
    if tree.is_tip[Int(node)]
        tip = tree.node_labels[Int(node)]
        return tip, tip
    end
    length(tree.children[Int(node)]) >= 2 || throw(ArgumentError("internal node $node has fewer than two child lineages"))
    left = _phyloref_first_tip_label(tree, map, tree.children[Int(node)][1])
    right = _phyloref_first_tip_label(tree, map, tree.children[Int(node)][2])
    return left <= right ? (left, right) : (right, left)
end

function phylo_branch_anchor(tree::CompactTree, edge::Integer; map::PhyloRef = build_phyloref(tree))
    1 <= edge <= tree.nedges || throw(ArgumentError("edge $edge is outside 1:$(tree.nedges)"))
    return phylo_node_anchor(tree, tree.child_of_edge[Int(edge)]; map = map)
end

function phylo_edge_signature(tree::CompactTree, edge::Integer; sep::AbstractString = "|", map::PhyloRef = build_phyloref(tree))
    1 <= edge <= tree.nedges || throw(ArgumentError("edge $edge is outside 1:$(tree.nedges)"))
    labels = [tree.tip_labels[Int(i)] for i in map.descendant_tip_positions_by_edge[Int(edge)]]
    sort!(labels)
    return join(labels, sep)
end

function phylo_edge_signatures(tree::CompactTree, edges::AbstractVector{<:Integer}; sep::AbstractString = "|", map::PhyloRef = build_phyloref(tree))
    return [phylo_edge_signature(tree, edge; sep = sep, map = map) for edge in edges]
end

function phylo_edges_from_signatures(tree::CompactTree, signatures::AbstractVector{<:AbstractString}; sep::AbstractString = "|", map::PhyloRef = build_phyloref(tree))
    index = Dict{String, Int}()
    for edge in 1:tree.nedges
        signature = phylo_edge_signature(tree, edge; sep = sep, map = map)
        haskey(index, signature) && throw(ArgumentError("edge signature is not unique: $signature"))
        index[signature] = edge
    end
    return [begin
        labels = split(String(signature), sep; keepempty = false)
        isempty(labels) && throw(ArgumentError("empty edge signature"))
        normalized = join(sort!(labels), sep)
        edge = get(index, normalized, 0)
        edge == 0 && throw(ArgumentError("no edge in tree has descendant-tip signature: $signature"))
        edge
    end for signature in signatures]
end

function R_node_id(tree::CompactTree, node::Integer; map::PhyloRef = build_phyloref(tree))
    1 <= node <= tree.nnodes || throw(ArgumentError("node $node is outside 1:$(tree.nnodes)"))
    return map.ape_node_ids[Int(node)]
end

function evotraits_node_id(tree::CompactTree, R_node_id_value::Integer; map::PhyloRef = build_phyloref(tree))
    1 <= R_node_id_value <= tree.nnodes || throw(ArgumentError("R node id $R_node_id_value is outside 1:$(tree.nnodes)"))
    node = map.evotraits_node_ids_from_ape[Int(R_node_id_value)]
    node != 0 || throw(ArgumentError("R node id $R_node_id_value has no EvoTraits mapping"))
    return node
end

function R_edge_id_cladewise(tree::CompactTree, edge::Integer; map::PhyloRef = build_phyloref(tree))
    1 <= edge <= tree.nedges || throw(ArgumentError("edge $edge is outside 1:$(tree.nedges)"))
    return map.ape_cladewise_edge_ranks[Int(edge)]
end

function R_edge_id_postorder(tree::CompactTree, edge::Integer; map::PhyloRef = build_phyloref(tree))
    1 <= edge <= tree.nedges || throw(ArgumentError("edge $edge is outside 1:$(tree.nedges)"))
    return map.ape_postorder_edge_ranks[Int(edge)]
end

function evotraits_edge_id_from_R_cladewise(tree::CompactTree, R_edge_id_value::Integer; map::PhyloRef = build_phyloref(tree))
    1 <= R_edge_id_value <= tree.nedges || throw(ArgumentError("R cladewise edge id $R_edge_id_value is outside 1:$(tree.nedges)"))
    return map.evotraits_edge_ids_from_ape_cladewise[Int(R_edge_id_value)]
end

function evotraits_edge_id_from_R_postorder(tree::CompactTree, R_edge_id_value::Integer; map::PhyloRef = build_phyloref(tree))
    1 <= R_edge_id_value <= tree.nedges || throw(ArgumentError("R postorder edge id $R_edge_id_value is outside 1:$(tree.nedges)"))
    return map.evotraits_edge_ids_from_ape_postorder[Int(R_edge_id_value)]
end

function R_edge_matrix(tree::CompactTree; order::Symbol = :cladewise, map::PhyloRef = build_phyloref(tree))
    ranks = order === :cladewise ? map.ape_cladewise_edge_ranks : order === :postorder ? map.ape_postorder_edge_ranks : throw(ArgumentError("Unsupported order=$order; expected :cladewise or :postorder"))
    edge_ids = sortperm(ranks)
    matrix = Matrix{Int}(undef, tree.nedges, 2)
    @inbounds for (i, edge) in enumerate(edge_ids)
        matrix[i, 1] = map.ape_node_ids[Int(tree.parent_of_edge[edge])]
        matrix[i, 2] = map.ape_node_ids[Int(tree.child_of_edge[edge])]
    end
    return matrix
end
