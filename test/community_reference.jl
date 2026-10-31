# Deliberately simple oracle. Descendant sets are built bottom-up, independently
# of production label matching and the ancestor walk used by the sparse kernels.
function community_reference(community, tree; tip_labels = tree.tip_labels)
    columns = Dict(label => col for (col, label) in enumerate(tip_labels))
    descendants = [Int[] for _ in 1:tree.nnodes]
    for node in tree.postorder
        if tree.is_tip[node]
            push!(descendants[node], columns[tree.node_labels[node]])
        else
            for child in tree.children[node]
                append!(descendants[node], descendants[child])
            end
        end
    end
    n = size(community, 1)
    counts = zeros(Int, n, tree.nedges)
    for edge in 1:tree.nedges, site in 1:n
        counts[site, edge] = count(col -> community[site, col] > 0,
            descendants[tree.child_of_edge[edge]])
    end
    active = counts .> 0
    occupancy = vec(sum(active; dims = 1))
    species_occupancy = vec(sum(community .> 0; dims = 1))
    richness = vec(sum(community .> 0; dims = 2))
    pd_values = zeros(n)
    pd_no_root = zeros(n)
    we, pe, strict, equal_pd = (zeros(n) for _ in 1:4)
    for site in 1:n
        for col in axes(community, 2)
            community[site, col] > 0 && (we[site] += 1 / species_occupancy[col])
        end
        for edge in 1:tree.nedges
            active[site, edge] || continue
            len = tree.edge_length[edge]
            pd_values[site] += len
            pe[site] += len / occupancy[edge]
            occupancy[edge] == 1 && (strict[site] += len)
            len > 0 && (equal_pd[site] += 1)
            counts[site, edge] < richness[site] && (pd_no_root[site] += len)
        end
        richness[site] == 1 && (pd_no_root[site] = NaN)
    end
    total = sum(tree.edge_length)
    positive = count(>(0), tree.edge_length)
    normalized = total > 0 ? pd_values ./ total : fill(NaN, n)
    equal_normalized = positive > 0 ? equal_pd ./ positive : fill(NaN, n)
    rpd = normalized ./ equal_normalized
    function beta(family)
        total, turnover, nestedness = (fill(NaN, n, n) for _ in 1:3)
        for i in 1:n, j in 1:n
            a = sum(tree.edge_length[active[i, :] .& active[j, :]])
            b = sum(tree.edge_length[active[i, :] .& .!active[j, :]])
            c = sum(tree.edge_length[.!active[i, :] .& active[j, :]])
            m = min(b, c)
            denom = family == :sorensen ? 2a + b + c : a + b + c
            turnover_denom = family == :sorensen ? a + m : a + 2m
            total[i, j] = denom > 0 ? (b + c) / denom : NaN
            turnover[i, j] = turnover_denom > 0 ?
                (family == :sorensen ? m : 2m) / turnover_denom : NaN
            nestedness[i, j] = total[i, j] - turnover[i, j]
        end
        return (; total, turnover, nestedness)
    end
    return (; active, occupancy, species_occupancy, richness, pd = pd_values,
        pd_no_root, we, pe, strict, normalized, equal_normalized, rpd, beta)
end

function random_community_tree(rng, ntips; zero_edges = true)
    len() = zero_edges && rand(rng) < 0.15 ? 0.0 : 0.1 + 3rand(rng)
    parts = ["s$i:$(len())" for i in 1:ntips]
    while length(parts) > 1
        shuffle!(rng, parts)
        width = rand(rng, 2:min(3, length(parts)))
        children = [pop!(parts) for _ in 1:width]
        push!(parts, "($(join(children, ','))):$(len())")
    end
    return only(parts) * ";"
end
