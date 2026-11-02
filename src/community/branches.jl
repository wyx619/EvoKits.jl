# Tree-dependent branch traversal is separate from community input handling.
function _community_roots(rows, tree::CompactTree)
    roots = fill(tree.root, length(rows))
    any(tips -> length(tips) > 1, rows) || return roots
    index = _build_tree_query_index(tree)
    for (site, tips) in enumerate(rows)
        isempty(tips) && continue
        ancestor = tree.tip_ids[first(tips)]
        for tip in Iterators.drop(tips, 1)
            ancestor = _tree_lca(index, ancestor, tree.tip_ids[tip])
        end
        roots[site] = ancestor
    end
    return roots
end

function _community_edge_index(tree::CompactTree)
    edge_of_child = zeros(Int32, tree.nnodes)
    for edge in 1:tree.nedges
        edge_of_child[tree.child_of_edge[edge]] = Int32(edge)
    end
    return edge_of_child
end

# A branch is visited once per site. Stop at the first visited ancestor, so
# neither a dense site-by-edge matrix nor all tip-to-root paths are needed.
function _walk_community_branches(f, rows, tree::CompactTree;
    roots::Union{Nothing, Vector{Int32}} = nothing,
    edge_of_child::Vector{Int32} = _community_edge_index(tree),
)
    seen = zeros(Int, tree.nnodes)
    @inbounds for (site, tips) in enumerate(rows)
        root = roots === nothing ? tree.root : roots[site]
        for tip in tips
            node = tree.tip_ids[tip]
            while node != root && seen[node] != site
                seen[node] = site
                edge = edge_of_child[node]
                f(site, edge)
                node = tree.parent_of_node[node]
            end
        end
    end
    return nothing
end

function _community_branch_sets(rows, tree::CompactTree)
    branches = [Int32[] for _ in eachindex(rows)]
    _walk_community_branches(rows, tree) do site, edge
        push!(branches[site], edge)
    end
    foreach(sort!, branches)
    return branches
end
