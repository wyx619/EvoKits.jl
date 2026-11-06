Base.@kwdef struct ExpectedPDResult
    richness::Vector{Int} = Int[]
    expected_pd::Vector{Float64} = Float64[]
    variance_pd::Vector{Float64} = Float64[]
    upper_bound::Bool = true
end

function edge_abundance_distribution(tree::CompactTree)
    counts = _tools_descendant_tip_counts(tree)
    return (num_children = Int[counts[Int(tree.child_of_edge[e])] for e in 1:tree.nedges],
        edge_length = copy(tree.edge_length))
end

ead(tree::CompactTree) = edge_abundance_distribution(tree)

@inline function _expected_pd_variance_term(length, absent, below, log_absent, factor)
    value = (length * absent) * below
    if absent >= floatmin(Float64) && value >= floatmin(Float64)
        return factor * value
    end
    (iszero(length) || iszero(below) || log_absent < -2200) && return 0.0
    # An absence probability can underflow even though branch-length products
    # make the final variance representable. Fall back only for these tails.
    return exp(log(length) + log(below) + log_absent + log(factor))
end

function _expected_pd_richness(richness, n::Int)
    richness === nothing && return collect(1:n)
    values = if richness isa Integer
        (richness,)
    elseif richness isa AbstractVector
        richness
    else
        throw(ArgumentError("richness must be an integer or a vector/range of integers"))
    end
    selected = Vector{Int}(undef, length(values))
    for (i, value) in enumerate(values)
        (value isa Integer && !(value isa Bool)) ||
            throw(ArgumentError("richness values must be integers, excluding Bool"))
        1 <= value <= n || throw(ArgumentError("richness values must lie in 1:$n"))
        selected[i] = Int(value)
    end
    return selected
end

function _expected_pd_vectors(tree::CompactTree, richness::Vector{Int}; upper_bound::Bool = true)
    _validate_community_branch_lengths(tree)
    n = tree.ntips
    expected, variance = zeros(length(richness)), zeros(length(richness))
    isempty(richness) && return expected, variance
    total_length = sum(tree.edge_length)
    if all(==(n), richness)
        fill!(expected, total_length)
        return expected, variance
    end
    descendants = _tools_descendant_tip_counts(tree)
    lengths = zeros(tree.nnodes)
    for edge in 1:tree.nedges
        lengths[tree.child_of_edge[edge]] = tree.edge_length[edge]
    end
    subtree = zeros(tree.nnodes)
    if upper_bound
        for node in tree.postorder
            parent = tree.parent_of_node[node]
            parent == 0 && continue
            subtree[parent] += subtree[node] + lengths[node]
        end
    end
    # Reuse probabilities for repeated descendant counts. No pair array or
    # LCA index is needed; all preparation and working arrays are linear.
    classes = unique(Int[descendants[node] for node in tree.child_of_edge])
    absent, present = zeros(n), zeros(n)
    for (i, count_tips) in enumerate(richness)
        if count_tips == n
            expected[i] = total_length
            continue
        end
        log_absence = log1p(-count_tips / n)
        for count in classes
            exponent = count * log_absence
            absent[count], present[count] = exp(exponent), -expm1(exponent)
        end
        upper_bound || fill!(subtree, 0.0)
        mean, var = 0.0, 0.0
        for node in tree.postorder
            parent = tree.parent_of_node[node]
            parent == 0 && continue
            count, length = Int(descendants[node]), lengths[node]
            a, q = absent[count], present[count]
            contribution = length * q
            mean += contribution
            # Multiply probabilities before squaring large branch lengths.
            exponent = count * log_absence
            var += _expected_pd_variance_term(length, a, contribution, exponent, 1.0)
            below = upper_bound ? subtree[node] * q : subtree[node]
            var += _expected_pd_variance_term(length, a, below, exponent, 2.0)
            upper_bound || (subtree[parent] += subtree[node] + contribution)
        end
        expected[i], variance[i] = mean, var
    end
    return expected, variance
end

"""
    expected_pd(tree; richness=nothing, upper_bound=true)
    pd_variance(tree; richness=nothing, upper_bound=true)

Expected rooted PD and its variance when each tip independently occurs with
probability richness/ntips. By default, return the full curve for richness
1:ntips. Supply an integer or a vector/range of integers in 1:ntips to compute
only those points, preserving request order and duplicates. Bool and
noninteger values are rejected. An empty vector returns empty output arrays;
a scalar still returns one-element arrays. Richness is the expected number
of tips, not a fixed-size sample without replacement.

`upper_bound=false` gives the exact Bernoulli variance; `true` keeps picante's
upper bound by using ancestor instead of descendant occurrence probabilities
in covariance terms. All-tip inclusion has variance zero.

Subtree sums replace ancestor/descendant edge-pair enumeration: O(nodes+ntips)
preparation plus O(m*nodes) work for m requested points, with
O(nodes+ntips+m) storage and no VCV, distance or edge-pair matrix. Branch
lengths must be finite/nonnegative with a finite sum, even for an empty
request. Inputs are not mutated. Returns `ExpectedPDResult` with owned arrays.
"""
function expected_pd(tree::CompactTree; richness = nothing, upper_bound::Bool = true)
    selected = _expected_pd_richness(richness, tree.ntips)
    expected, variance = _expected_pd_vectors(tree, selected; upper_bound)
    return ExpectedPDResult(richness = selected, expected_pd = expected,
        variance_pd = variance, upper_bound = upper_bound)
end

pd_variance(tree::CompactTree; richness = nothing, upper_bound::Bool = true) =
    expected_pd(tree; richness, upper_bound)
