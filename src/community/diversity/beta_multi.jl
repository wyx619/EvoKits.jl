"""
    BetaMultiResult

Multiple-site taxonomic beta decomposition from [`beta_multi`](@ref).
`total`, `turnover` and `nestedness` are scalars for the selected collection.
Owned `sample_labels`, `sample_indices` and `richness` follow the requested
site order; indices refer to the original input rows. `pooled_richness` counts
species present in at least one selected site. `index_family` records the family.
"""
struct BetaMultiResult
    sample_labels::Vector{String}
    sample_indices::Vector{Int}
    richness::Vector{Int}
    pooled_richness::Int
    total::Float64
    turnover::Float64
    nestedness::Float64
    index_family::Symbol
end

function _beta_multi_samples(samples, labels)
    samples === nothing && return collect(eachindex(labels))
    samples isa Union{AbstractVector, Tuple} ||
        throw(ArgumentError("samples must be a vector or tuple of row indices or string labels"))
    indices = if all(x -> x isa Integer && !(x isa Bool), samples)
        all(i -> 1 <= i <= length(labels), samples) ||
            throw(ArgumentError("samples contains a row index outside the input"))
        Int[Int(i) for i in samples]
    elseif all(x -> x isa AbstractString, samples)
        positions = Dict(label => i for (i, label) in enumerate(labels))
        all(label -> haskey(positions, label), samples) ||
            throw(ArgumentError("samples contains an unknown sample label"))
        Int[positions[label] for label in samples]
    else
        throw(ArgumentError("samples must contain only integer row indices or only string labels; Bool masks are not supported"))
    end
    allunique(indices) || throw(ArgumentError("samples must not contain duplicate sites"))
    return indices
end

function _beta_multi_counts(community::AbstractMatrix{<:Real}, indices)
    richness, occupancy = zeros(Int, length(indices)), zeros(Int, size(community, 2))
    for col in axes(community, 2), (site, row) in enumerate(indices)
        if community[row, col] > 0
            richness[site] += 1
            occupancy[col] += 1
        end
    end
    return richness, occupancy
end

function _beta_multi_counts(community::SparseMatrixCSC{<:Real}, indices)
    positions = zeros(Int, size(community, 1))
    for (site, row) in enumerate(indices)
        positions[row] = site
    end
    richness, occupancy = zeros(Int, length(indices)), zeros(Int, size(community, 2))
    rows, values = rowvals(community), nonzeros(community)
    for col in axes(community, 2), entry in nzrange(community, col)
        site = positions[rows[entry]]
        if site > 0 && values[entry] > 0
            richness[site] += 1
            occupancy[col] += 1
        end
    end
    return richness, occupancy
end

# For occupancy q, sum of pairwise shared counts is q(q-1)/2. For sorted
# richness r[k], pairwise minima/maxima sum to sum((n-k)r[k]) and
# sum((k-1)r[k]), respectively.
# Subtract shared counts to obtain min/max unique counts, without visiting
# pairs. Int128 keeps these cancellations exact, including totals above 2^53.
function _beta_multi_statistics(richness, occupancy)
    shared_excess, shared_pairs, pooled = Int128(0), Int128(0), 0
    for q in occupancy
        if q > 0
            shared_excess += q - 1
            shared_pairs += Int128(q) * (q - 1) ÷ 2
            pooled += 1
        end
    end
    minimum_unique = maximum_unique = -shared_pairs
    n = length(richness)
    for (k, r) in enumerate(sort(richness))
        minimum_unique += Int128(n - k) * r
        maximum_unique += Int128(k - 1) * r
    end
    return (; shared_excess, minimum_unique, maximum_unique, pooled)
end

function _beta_multi_components(a, minimum_unique, maximum_unique, family)
    unique = minimum_unique + maximum_unique
    total_denom = family === :sorensen ? 2a + unique : a + unique
    turnover_num = family === :sorensen ? minimum_unique : 2minimum_unique
    turnover_denom = a + turnover_num
    total = total_denom > 0 ? unique / total_denom : NaN
    turnover = turnover_denom > 0 ? turnover_num / turnover_denom : NaN
    # Use the defining product, retaining small nestedness when total and
    # turnover round to the same Float64. Exact integer differences come first.
    nestedness = total_denom > 0 && turnover_denom > 0 ?
        (a / turnover_denom) * ((maximum_unique - minimum_unique) / total_denom) : NaN
    return total, turnover, nestedness
end

"""
    beta_multi(community; index_family=:sorensen, samples=nothing, sample_labels=nothing)

Decompose presence-based taxonomic beta diversity of an entire collection
into scalar total, turnover and nestedness components. Supports `:sorensen`
and `:jaccard`, dense/CSC matrices and `CommunityData`; positive abundance
denotes presence. This is a multiple-site index, not a mean of pairwise beta.

`samples` selects distinct rows by integer indices or string labels in a vector
or tuple. Omit it to select every row. At least two selected sites are required;
Bool masks, repeated sites and mixed index/label selectors are errors. Integer
sample identifiers become strings in `sample_labels`; select them by string
when using labels. Labels must be unique and nonempty. Results retain selection
order, copy labels and indices, and count pooled species only within that subset.
The complete input must contain finite non-negative real values, including
unselected rows; stored zeros are absent. Inputs are unchanged.

Let A = sum of selected richness minus pooled richness, and M/X be the sums
of minimum/maximum unique species counts over unordered site pairs. Sørensen
total is (M+X)/(2A+M+X) and turnover is M/(A+M). Jaccard total is
(M+X)/(A+M+X) and turnover is 2M/(A+2M). Nestedness uses the equivalent
defining product to avoid subtracting nearly equal components.

Empty sites are retained. All-empty input gives three NaNs. If no species is
shared among sites, total is one when any species is present; turnover is one
when M>0 and otherwise NaN, with nestedness zero or NaN respectively.
Other zero denominators yield NaN for the affected components.

Occupancy and sorted richness replace pair enumeration: CSC work is
O(nnz + p + n + m log m), with O(n+p) auxiliary storage for n input sites,
p species and m selected sites. No dense community or site-pair matrix is built.
"""
function beta_multi(community::AbstractMatrix{<:Real};
    index_family::Symbol = :sorensen, samples = nothing, sample_labels = nothing,
)
    _validate_community_beta_family(index_family)
    _validate_community_values(community)
    labels = _validated_labels(sample_labels, size(community, 1), "sample")
    indices = _beta_multi_samples(samples, labels)
    length(indices) >= 2 || throw(ArgumentError("beta_multi requires at least two selected sites"))
    richness, occupancy = _beta_multi_counts(community, indices)
    stats = _beta_multi_statistics(richness, occupancy)
    total, turnover, nestedness = _beta_multi_components(stats.shared_excess,
        stats.minimum_unique, stats.maximum_unique, index_family)
    return BetaMultiResult(labels[indices], indices, richness, stats.pooled,
        total, turnover, nestedness, index_family)
end

beta_multi(community::AbstractMatrix; kwargs...) =
    throw(ArgumentError("community matrix must contain real numbers without missing values"))
