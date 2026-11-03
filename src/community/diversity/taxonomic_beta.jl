"""
    BetaDissResult

Taxonomic beta decomposition returned by [`beta_diss`](@ref). `total`,
`turnover`, and `nestedness` are symmetric matrices, or vectors when `pairs`
is supplied. `sample_labels` describe the complete input row order; `pairs`
contains the requested row indices, or `nothing` for full matrices.
"""
struct BetaDissResult{T<:AbstractArray{Float64}}
    sample_labels::Vector{String}
    total::T
    turnover::T
    nestedness::T
    index_family::Symbol
    pairs::Union{Nothing, Vector{Tuple{Int, Int}}}
end

"""
    beta_diss(community; index_family=:sorensen, pairs=nothing,
              sample_labels=nothing)

Partition presence-based taxonomic beta diversity into turnover and nestedness
without a phylogenetic tree. Rows are sites and columns are species. Values
must be finite and non-negative; every positive value denotes presence,
irrespective of abundance. Dense matrices and `SparseMatrixCSC` are supported
without constructing a dense site-by-species matrix.

For shared species count `a` and unique counts `b`, `c`, Sørensen total is
`(b+c)/(2a+b+c)` and turnover is `min(b,c)/(a+min(b,c))`. Jaccard total is
`(b+c)/(a+b+c)` and turnover is `2min(b,c)/(a+2min(b,c))`. Nestedness is
total minus turnover. Supports `:sorensen` and `:jaccard`.

Returns a `BetaDissResult` with three symmetric matrices in input row order.
Supply `pairs` as an iterable of `(row_i, row_j)` to return three vectors in
the supplied order, retaining repeated, reversed, and self pairs. Selected
pairs avoid constructing full site-by-site matrices. Labels are copied.

Zero denominators return `NaN`. An empty site compared with a nonempty site
has total one and undefined turnover/nestedness. Two empty sites have `NaN`
for all components, including the diagonal. Nonempty diagonal entries are zero.
"""
function beta_diss(community::AbstractMatrix{<:Real};
    index_family::Symbol = :sorensen,
    pairs = nothing,
    sample_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
)
    _validate_community_beta_family(index_family)
    n, nspecies = size(community)
    labels = _tools_sample_labels(n, sample_labels)
    selected = _community_pairs(pairs, n)
    rows = _community_presence_rows(community, axes(community, 2))
    total, turnover, nestedness = _taxonomic_beta(rows, nspecies, index_family, selected)
    return BetaDissResult(labels, total, turnover, nestedness, index_family, selected)
end

# Species rows are sorted and unique by construction. Counts remain exact;
# unlike branch lengths, subtracting shared counts cannot lose small weights.
function _community_species_components(left, right)
    a, i, j = 0, 1, 1
    @inbounds while i <= length(left) && j <= length(right)
        if left[i] == right[j]
            a += 1
            i += 1
            j += 1
        elseif left[i] < right[j]
            i += 1
        else
            j += 1
        end
    end
    return a, length(left) - a, length(right) - a
end

function _taxonomic_beta(rows, nspecies, index_family, selected)
    n = length(rows)
    if selected === nothing && n >= 128
        incidence = _community_incidence(rows, nspecies)
        shared = Matrix(transpose(incidence) * incidence)
        richness = length.(rows)
        return _community_beta(n, index_family, selected; total = shared) do i, j
            a = shared[i, j]
            (a, richness[i] - a, richness[j] - a)
        end
    end
    return _community_beta(n, index_family, selected) do i, j
        _community_species_components(rows[i], rows[j])
    end
end
