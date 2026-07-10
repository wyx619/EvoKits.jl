Base.@kwdef struct PhyloRef
    tree::CompactTree
    ape_node_ids::Vector{Int} = Int[]
    evotraits_node_ids_from_ape::Vector{Int} = Int[]
    ape_cladewise_edge_ranks::Vector{Int} = Int[]
    ape_postorder_edge_ranks::Vector{Int} = Int[]
    evotraits_edge_ids_from_ape_cladewise::Vector{Int} = Int[]
    evotraits_edge_ids_from_ape_postorder::Vector{Int} = Int[]
    tip_positions_by_node::Vector{Int} = Int[]
    descendant_tip_positions_by_node::Vector{Vector{Int32}} = Vector{Vector{Int32}}()
    descendant_tip_positions_by_edge::Vector{Vector{Int32}} = Vector{Vector{Int32}}()
end

function build_phyloref(tree::CompactTree)
    ape_node_ids = fill(0, tree.nnodes)
    evotraits_node_ids_from_ape = fill(0, tree.nnodes)
    tip_positions_by_node = fill(0, tree.nnodes)

    @inbounds for (tip_pos, node) in enumerate(tree.tip_ids)
        ape_node_ids[Int(node)] = tip_pos
        evotraits_node_ids_from_ape[tip_pos] = Int(node)
        tip_positions_by_node[Int(node)] = tip_pos
    end

    preorder_internal = Int32[node for node in tree.preorder if !tree.is_tip[Int(node)]]
    @inbounds for (internal_rank, node) in enumerate(preorder_internal)
        ape_id = tree.ntips + internal_rank
        ape_node_ids[Int(node)] = ape_id
        evotraits_node_ids_from_ape[ape_id] = Int(node)
    end

    descendant_tip_positions_by_node = [Int32[] for _ in 1:tree.nnodes]
    @inbounds for node in tree.postorder
        inode = Int(node)
        if tree.is_tip[inode]
            push!(descendant_tip_positions_by_node[inode], Int32(tip_positions_by_node[inode]))
        else
            tips = descendant_tip_positions_by_node[inode]
            for child in tree.children[inode]
                append!(tips, descendant_tip_positions_by_node[Int(child)])
            end
            sort!(tips)
        end
    end

    descendant_tip_positions_by_edge = Vector{Vector{Int32}}(undef, tree.nedges)
    @inbounds for edge in 1:tree.nedges
        descendant_tip_positions_by_edge[edge] = copy(descendant_tip_positions_by_node[Int(tree.child_of_edge[edge])])
    end

    cladewise_edge_ids = Int[]
    function append_cladewise_edges!(node::Int)
        for edge in Int(tree.first_child_edge[node]):Int(tree.last_child_edge[node])
            push!(cladewise_edge_ids, edge)
            child = Int(tree.child_of_edge[edge])
            tree.is_tip[child] || append_cladewise_edges!(child)
        end
    end
    append_cladewise_edges!(Int(tree.root))

    postorder_edge_ids = Int[]
    function append_postorder_edges!(node::Int)
        for edge in Int(tree.last_child_edge[node]):-1:Int(tree.first_child_edge[node])
            child = Int(tree.child_of_edge[edge])
            tree.is_tip[child] || append_postorder_edges!(child)
        end
        append!(postorder_edge_ids, Int(tree.first_child_edge[node]):Int(tree.last_child_edge[node]))
    end
    append_postorder_edges!(Int(tree.root))

    ape_cladewise_edge_ranks = fill(0, tree.nedges)
    ape_postorder_edge_ranks = fill(0, tree.nedges)
    evotraits_edge_ids_from_ape_cladewise = Vector{Int}(undef, tree.nedges)
    evotraits_edge_ids_from_ape_postorder = Vector{Int}(undef, tree.nedges)
    @inbounds for (rank, edge) in enumerate(cladewise_edge_ids)
        ape_cladewise_edge_ranks[edge] = rank
        evotraits_edge_ids_from_ape_cladewise[rank] = edge
    end
    @inbounds for (rank, edge) in enumerate(postorder_edge_ids)
        ape_postorder_edge_ranks[edge] = rank
        evotraits_edge_ids_from_ape_postorder[rank] = edge
    end

    return PhyloRef(
        tree = tree,
        ape_node_ids = ape_node_ids,
        evotraits_node_ids_from_ape = evotraits_node_ids_from_ape,
        ape_cladewise_edge_ranks = ape_cladewise_edge_ranks,
        ape_postorder_edge_ranks = ape_postorder_edge_ranks,
        evotraits_edge_ids_from_ape_cladewise = evotraits_edge_ids_from_ape_cladewise,
        evotraits_edge_ids_from_ape_postorder = evotraits_edge_ids_from_ape_postorder,
        tip_positions_by_node = tip_positions_by_node,
        descendant_tip_positions_by_node = descendant_tip_positions_by_node,
        descendant_tip_positions_by_edge = descendant_tip_positions_by_edge,
    )
end
