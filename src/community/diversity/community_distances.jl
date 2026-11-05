# Store only observed tips and weights. Scale each row before normalizing so
# finite abundances do not overflow in row sums or pairwise outer products.
@inline function _community_distance_ratio(value, scale)
    numerator, denominator = Float64(value), Float64(scale)
    if isfinite(numerator) && isfinite(denominator) && denominator > 0 &&
        (numerator > 0 || iszero(value))
        return numerator / denominator
    end
    # Large integers/BigFloats may not fit Float64 even when their ratio does.
    return Float64(BigFloat(value) / BigFloat(scale))
end
_community_distance_ratio(value::Rational, scale::Rational) =
    Float64(BigFloat(value) / BigFloat(scale))

function _community_distance_rows(community, tree, tip_labels, weighted;
    active_sites = nothing, need_lengths = true)
    _validate_community_values(community)
    _validate_community_branch_lengths(tree)
    size(community, 2) == tree.ntips ||
        throw(ArgumentError("community matrix must have one column per tree tip"))
    perm = _tools_tip_permutation(tree, tip_labels)
    n = size(community, 1)
    nodes = [Int32[] for _ in 1:n]
    T = weighted ? eltype(community) : Float64
    raw = [T[] for _ in 1:n]
    _community_each_positive(community) do row, col, value
        active_sites !== nothing && !active_sites[row] && return
        push!(nodes[row], tree.tip_ids[perm[col]])
        push!(raw[row], weighted ? value : 1.0)
    end
    weights = [Float64[] for _ in 1:n]
    scales = zeros(T, n)
    totals = zeros(n)
    for row in 1:n
        isempty(raw[row]) && continue
        scales[row] = maximum(raw[row])
        weights[row] = [_community_distance_ratio(x, scales[row]) for x in raw[row]]
        totals[row] = sum(weights[row])
        weights[row] ./= totals[row]
    end
    lengths = need_lengths ? _comdist_lengths(tree) : Float64[]
    return (; nodes, weights, raw, scales, totals, lengths)
end

function _comdist_lengths(tree)
    lengths = zeros(tree.nnodes)
    for edge in 1:tree.nedges
        lengths[tree.child_of_edge[edge]] = tree.edge_length[edge]
    end
    return lengths
end

# The distance at a node is the target-weighted sum below it plus the sum
# outside its subtree. Prefix/suffix sibling sums avoid subtracting nearly
# equal masses or long shared root paths.
struct _ComdistWorkspace
    mass::Vector{Float64}
    down::Vector{Float64}
    outside_mass::Vector{Float64}
    outside_distance::Vector{Float64}
    suffix_mass::Vector{Float64}
    suffix_distance::Vector{Float64}
end

_ComdistWorkspace(n) = _ComdistWorkspace((zeros(n) for _ in 1:6)...)

function _comdist_target!(workspace, tree, nodes, weights, lengths)
    mass, down = workspace.mass, workspace.down
    fill!(mass, 0)
    fill!(down, 0)
    for (node, weight) in zip(nodes, weights)
        mass[node] = weight
    end
    for node in tree.postorder
        parent = tree.parent_of_node[node]
        parent == 0 && continue
        mass[parent] += mass[node]
        down[parent] += down[node] + lengths[node] * mass[node]
    end
    outside_mass, outside_distance = workspace.outside_mass, workspace.outside_distance
    outside_mass[tree.root] = outside_distance[tree.root] = 0.0
    for parent in tree.preorder
        after_mass, after_distance = 0.0, 0.0
        for child in Iterators.reverse(tree.children[parent])
            workspace.suffix_mass[child] = after_mass
            workspace.suffix_distance[child] = after_distance
            after_mass += mass[child]
            after_distance += down[child] + lengths[child] * mass[child]
        end
        before_mass, before_distance = 0.0, 0.0
        for child in tree.children[parent]
            outside_mass[child] = outside_mass[parent] + before_mass + workspace.suffix_mass[child]
            outside_distance[child] = outside_distance[parent] + before_distance +
                workspace.suffix_distance[child] + lengths[child] * outside_mass[child]
            before_mass += mass[child]
            before_distance += down[child] + lengths[child] * mass[child]
        end
    end
    return nothing
end

"""
    comdist(community, tree; abundance_weighted=false, pairs=nothing, tip_labels=nothing)
    comdist(data::CommunityData, tree; abundance_weighted=false, pairs=nothing)

Mean pairwise phylogenetic distance between communities. Rows are sites and
columns are tree tips; `tip_labels` describes matrix column order. Positive
values denote presence, or supply weights with `abundance_weighted=true`.
Returns a symmetric Float64 matrix in site order by default. With `pairs`,
returns a Float64 vector in the supplied order of integer site-index pairs,
retaining repeated, reversed and self pairs. Empty requests give an empty
vector. Input values/labels and all branch lengths are always validated.
Only requested endpoints are retained for selected-pair computation.
Diagonals include all
within-site ordered tip pairs, including zero self-distances; empty sites
give NaN. Branch lengths must be finite and non-negative with a finite sum.

CSC input is traversed directly. Sparse requests can compress the tree to
observed tips and connecting nodes, summing original branches on unary paths;
otherwise uses full-tree passes. Preparation is linear in tree size, without
a full tip-distance, dense community or site-pair weight matrix. Full output
remains quadratic in sites; selected output uses linear space in requests.
Inputs are not mutated.
"""
function comdist(community::AbstractMatrix{<:Real}, tree::CompactTree;
    abundance_weighted::Bool = false,
    pairs = nothing,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    selected = _community_pairs(pairs, size(community, 1))
    active_sites = selected === nothing ? nothing : falses(size(community, 1))
    if selected !== nothing
        for (i, j) in selected
            active_sites[i] = active_sites[j] = true
        end
    end
    rows = _community_distance_rows(community, tree, tip_labels, abundance_weighted;
        active_sites, need_lengths = false)
    return _comdist_adaptive(rows, tree, selected)
end

# RaoD reuses the prepared rows and workspace, then evaluates the pooled
# community with the same tree traversal. No input needs to be aligned twice.
function _comdist(rows, tree, weights = rows.weights, workspace = _ComdistWorkspace(tree.nnodes))
    n = length(rows.nodes)
    result = fill(NaN, n, n)
    for target in 1:n
        isempty(rows.nodes[target]) && continue
        _comdist_target!(workspace, tree, rows.nodes[target], weights[target], rows.lengths)
        for source in target:n
            isempty(rows.nodes[source]) && continue
            value = 0.0
            for (node, weight) in zip(rows.nodes[source], weights[source])
                value += weight * (workspace.down[node] + workspace.outside_distance[node])
            end
            result[source, target] = result[target, source] = value
        end
    end
    return result
end

# Keep the two closest DISTINCT target identities. This supports excluding
# the query species while retaining zero distances to other species.
struct _ComdistntWorkspace
    first_distance::Vector{Float64}
    second_distance::Vector{Float64}
    first_tip::Vector{Int32}
    second_tip::Vector{Int32}
end

_ComdistntWorkspace(n) = _ComdistntWorkspace(fill(Inf, n), fill(Inf, n), zeros(Int32, n), zeros(Int32, n))

@inline function _comdistnt_insert!(workspace, node, distance, tip)
    tip == 0 && return nothing
    d1, d2 = workspace.first_distance, workspace.second_distance
    t1, t2 = workspace.first_tip, workspace.second_tip
    if tip == t1[node]
        d1[node] = min(d1[node], distance)
    elseif tip == t2[node]
        d2[node] = min(d2[node], distance)
        if d2[node] < d1[node]
            d1[node], d2[node] = d2[node], d1[node]
            t1[node], t2[node] = t2[node], t1[node]
        end
    elseif distance < d1[node]
        d2[node], t2[node] = d1[node], t1[node]
        d1[node], t1[node] = distance, tip
    elseif distance < d2[node]
        d2[node], t2[node] = distance, tip
    end
    return nothing
end

function _comdistnt_target!(workspace, tree, nodes, lengths)
    fill!(workspace.first_distance, Inf)
    fill!(workspace.second_distance, Inf)
    fill!(workspace.first_tip, 0)
    fill!(workspace.second_tip, 0)
    for node in nodes
        workspace.first_distance[node] = 0.0
        workspace.first_tip[node] = node
    end
    for node in tree.postorder
        parent = tree.parent_of_node[node]
        parent == 0 && continue
        _comdistnt_insert!(workspace, parent,
            workspace.first_distance[node] + lengths[node], workspace.first_tip[node])
        _comdistnt_insert!(workspace, parent,
            workspace.second_distance[node] + lengths[node], workspace.second_tip[node])
    end
    for node in tree.preorder
        parent = tree.parent_of_node[node]
        parent == 0 && continue
        # Descendant targets reached via the parent have a longer backtracking
        # path; the child's existing downward candidates dominate those copies.
        _comdistnt_insert!(workspace, node,
            workspace.first_distance[parent] + lengths[node], workspace.first_tip[parent])
        _comdistnt_insert!(workspace, node,
            workspace.second_distance[parent] + lengths[node], workspace.second_tip[parent])
    end
    return nothing
end

function _comdistnt_directional(workspace, nodes, raw, exclude)
    # Scale after excluding ineligible sources: a dominant conspecific may
    # otherwise erase the tiny weights of all remaining eligible tips.
    total, scale = 0.0, zero(eltype(raw))
    for (node, weight) in zip(nodes, raw)
        distance = exclude && workspace.first_tip[node] == node ?
            workspace.second_distance[node] : workspace.first_distance[node]
        if isfinite(distance)
            if weight > scale
                total = total * _community_distance_ratio(scale, weight) + 1
                scale = weight
            else
                total += _community_distance_ratio(weight, scale)
            end
        end
    end
    total > 0 || return NaN
    value = 0.0
    for (node, weight) in zip(nodes, raw)
        distance = exclude && workspace.first_tip[node] == node ?
            workspace.second_distance[node] : workspace.first_distance[node]
        if isfinite(distance)
            value += distance * (_community_distance_ratio(weight, scale) / total)
        end
    end
    return value
end

function _comdistnt_symmetric(forward, reverse, rows, i, j, weighted)
    isfinite(forward) || return reverse
    isfinite(reverse) || return forward
    weighted || return forward / 2 + reverse / 2
    # Preserve the original site-total weighting between directions without
    # overflowing raw totals. Within directions only eligible query tips count.
    scale = max(rows.scales[i], rows.scales[j])
    wi = _community_distance_ratio(rows.scales[i], scale) * rows.totals[i]
    wj = _community_distance_ratio(rows.scales[j], scale) * rows.totals[j]
    total = wi + wj
    return forward * (wi / total) + reverse * (wj / total)
end

"""
    comdistnt(community, tree; abundance_weighted=false,
              exclude_conspecifics=false, tip_labels=nothing)
    comdistnt(data::CommunityData, tree; kwargs...)

Symmetric mean of the two directional nearest-taxon distances. Matrix rows
are sites; `tip_labels` describes the column order. Each direction averages
over present source species with a valid target. With `abundance_weighted=true`,
source abundances weight that mean, and original site totals weight the two
directions. Otherwise directions receive equal weight. If one direction has
no valid targets, use the other; if neither has targets, return NaN.

`exclude_conspecifics=true` excludes the same species identity, retaining
zero distances to different species. With exclusion, the diagonal is the
within-site nearest-other-species mean (NaN for a singleton). With the default
setting, nonempty diagonals are zero. Empty sites give NaN. Branch lengths
must be finite and non-negative with a finite sum.

Two tree passes per nonempty target propagate the nearest two distinct target
species, without building a full tip-distance matrix or dense community. CSC
input is traversed directly; output is an ordinary symmetric Float64 matrix
in site order. Neither the community nor the tree is mutated.
"""
function comdistnt(community::AbstractMatrix{<:Real}, tree::CompactTree;
    abundance_weighted::Bool = false,
    exclude_conspecifics::Bool = false,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    rows = _community_distance_rows(community, tree, tip_labels, abundance_weighted)
    n = length(rows.nodes)
    result = fill(NaN, n, n)
    workspace = _ComdistntWorkspace(tree.nnodes)
    for target in 1:n
        isempty(rows.nodes[target]) && continue
        _comdistnt_target!(workspace, tree, rows.nodes[target], rows.lengths)
        for source in 1:n
            forward = _comdistnt_directional(workspace,
                rows.nodes[source], rows.raw[source], exclude_conspecifics)
            if source >= target
                result[source, target] = forward
            else
                reverse = result[target, source]
                value = _comdistnt_symmetric(forward, reverse, rows, source, target, abundance_weighted)
                result[source, target] = result[target, source] = value
            end
        end
    end
    return result
end
