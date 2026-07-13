struct CompactTree
    ntips::Int
    nnodes::Int
    nedges::Int
    parent_of_edge::Vector{Int32}
    child_of_edge::Vector{Int32}
    edge_length::Vector{Float64}
    root::Int32
    parent_of_node::Vector{Int32}
    dist_from_root::Vector{Float64}
    is_tip::BitVector
    tip_ids::Vector{Int32}
    postorder::Vector{Int32}
    postorder_internal::Vector{Int32}
    preorder::Vector{Int32}
    children::Vector{Vector{Int32}}
    first_child_edge::Vector{Int32}
    last_child_edge::Vector{Int32}
    tipname_to_id::Dict{String, Int32}
    tip_labels::Vector{String}
    node_labels::Vector{String}
end
