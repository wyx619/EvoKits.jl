"""
    CommunityData(matrix; sample_labels=nothing, species_labels=nothing, storage=:preserve)
    CommunityData(table; sample=:sample, species=:species, value=:value,
                  sample_labels=nothing, species_labels=nothing, duplicates=:sum)
    CommunityData(rows, cols, values; sample_labels, species_labels, duplicates=:sum)

Labeled community snapshot with fields `matrix`, `sample_labels`, and
`species_labels`. Rows are sites and columns are species. Values must be finite
and non-negative. Abundances are preserved; presence-based metrics interpret
all positive values as presence. Treat the owned arrays as read-only.

Matrix input is copied, preserving its numeric and storage types by default.
Use `storage=:sparse` to construct CSC storage. Omitted matrix labels become
`sample_1`, ... and `species_1`, ...; provide actual species names for tree
alignment. Labels must be unique, nonempty strings or integer identifiers;
integer identifiers are converted to strings.

An `AbstractDataFrame` is interpreted as a long table. Select the site, species
and abundance columns with `sample`, `species`, and `value`; `value=nothing`
counts occurrences. Labels default to first-occurrence order, including rows
with zero abundance. Explicit complete label lists retain unobserved species
and empty sites, and must include all identifiers appearing in the table.

Long tables and indexed triplets build CSC storage directly. Triplets use
one-based integer row/column indices into the required complete label lists.
Duplicate coordinates are summed (`duplicates=:sum`) or rejected
(`duplicates=:error`), including duplicates with zero values. Bool records
are converted to integer counts. Invalid inputs and overflowing sums fail;
stored zeros are removed from constructed CSC matrices.

Use [`align_community`](@ref) for tree alignment, or pass the data directly to
community metrics and [`prepare_community`](@ref). Constructors copy their inputs;
results and prepared contexts may be modified independently of these inputs.
"""
struct CommunityData{M<:AbstractMatrix{<:Real}}
    matrix::M
    sample_labels::Vector{String}
    species_labels::Vector{String}

    function CommunityData(matrix::M, samples::Vector{String}, species::Vector{String},
        ::Val{:validated}) where {M<:AbstractMatrix{<:Real}}
        new{M}(matrix, samples, species)
    end
end

Base.size(data::CommunityData) = size(data.matrix)
Base.size(data::CommunityData, dim::Integer) = size(data.matrix, dim)

function _community_identifier(value, name)
    (value isa AbstractString || (value isa Integer && !(value isa Bool))) ||
        throw(ArgumentError("$name must be nonempty strings or integer identifiers"))
    label = String(string(value))
    isempty(strip(label)) && throw(ArgumentError("$name must not be empty"))
    return label
end

function _community_data_labels(labels, n, prefix)
    labels === nothing && return ["$(prefix)_$i" for i in 1:n]
    labels isa AbstractVector || throw(ArgumentError("$(prefix)_labels must be a vector"))
    length(labels) == n || throw(ArgumentError("$(prefix)_labels length must match the matrix dimension"))
    result = [_community_identifier(label, "$(prefix)_labels") for label in labels]
    length(unique(result)) == n || throw(ArgumentError("$(prefix)_labels must not contain duplicates"))
    return result
end

function CommunityData(matrix::AbstractMatrix{<:Real};
    sample_labels = nothing, species_labels = nothing, storage::Symbol = :preserve,
)
    _validate_community_values(matrix)
    storage in (:preserve, :sparse) || throw(ArgumentError("storage must be :preserve or :sparse"))
    samples = _community_data_labels(sample_labels, size(matrix, 1), "sample")
    species = _community_data_labels(species_labels, size(matrix, 2), "species")
    owned = matrix isa SparseMatrixCSC || storage === :preserve ? copy(matrix) : sparse(matrix)
    owned isa SparseMatrixCSC && dropzeros!(owned)
    return CommunityData(owned, samples, species, Val(:validated))
end

CommunityData(matrix::AbstractMatrix; kwargs...) =
    throw(ArgumentError("community matrix must contain real numbers without missing values"))

function _community_record_values(values)
    all(x -> x isa Real && isfinite(x) && x >= 0, values) ||
        throw(ArgumentError("community values must be finite, non-negative real numbers"))
    isempty(values) && return Float64[]
    # Concretize heterogeneous table columns without imposing Float64 on counts.
    T = foldl(promote_type, (x isa Bool ? Int : typeof(x) for x in values))
    result = T[x for x in values]
    all(isfinite, result) || throw(ArgumentError("community values overflow after numeric promotion"))
    return result
end

function _community_duplicate_sum(a::T, b::T) where {T}
    result = try
        T <: Integer && isbitstype(T) ? Base.Checked.checked_add(a, b) : a + b
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("duplicate abundance sum overflowed"))
    end
    isfinite(result) && result >= 0 || throw(ArgumentError("duplicate abundance sum must remain finite and non-negative"))
    return result
end

function _community_triplet_matrix(rows, cols, values, n, p, duplicates)
    duplicates in (:sum, :error) || throw(ArgumentError("duplicates must be :sum or :error"))
    length(rows) == length(cols) == length(values) ||
        throw(ArgumentError("rows, cols and values must have equal lengths"))
    all(i -> i isa Integer && !(i isa Bool) && 1 <= i <= n, rows) ||
        throw(ArgumentError("rows must contain valid one-based integer site indices"))
    all(j -> j isa Integer && !(j isa Bool) && 1 <= j <= p, cols) ||
        throw(ArgumentError("cols must contain valid one-based integer species indices"))
    numeric = _community_record_values(values)
    combine = duplicates === :sum ? _community_duplicate_sum :
        ((a, b) -> throw(ArgumentError("duplicate site/species records are not allowed")))
    matrix = sparse(Int.(rows), Int.(cols), numeric, n, p, combine)
    dropzeros!(matrix)
    return matrix
end

function CommunityData(rows::AbstractVector, cols::AbstractVector, values::AbstractVector;
    sample_labels::AbstractVector, species_labels::AbstractVector, duplicates::Symbol = :sum,
)
    Base.require_one_based_indexing(rows, cols, values)
    samples = _community_data_labels(sample_labels, length(sample_labels), "sample")
    species = _community_data_labels(species_labels, length(species_labels), "species")
    matrix = _community_triplet_matrix(rows, cols, values, length(samples), length(species), duplicates)
    return CommunityData(matrix, samples, species, Val(:validated))
end

function _community_table_column(table, column)
    (column isa Symbol || column isa AbstractString) && String(column) in names(table) ||
        throw(ArgumentError("long table is missing column: $column"))
    return table[!, column]
end

function _community_record_indices(values, labels, prefix)
    records = [_community_identifier(value, "$prefix identifiers") for value in values]
    domain = labels === nothing ? unique(records) :
        _community_data_labels(labels, length(labels), prefix)
    positions = Dict(label => i for (i, label) in enumerate(domain))
    indices = Vector{Int}(undef, length(records))
    for (i, label) in enumerate(records)
        index = get(positions, label, 0)
        index > 0 || throw(ArgumentError("$prefix identifier $label is absent from $(prefix)_labels"))
        indices[i] = index
    end
    return indices, domain
end

function CommunityData(table::AbstractDataFrame;
    sample = :sample, species = :species, value = :value,
    sample_labels = nothing, species_labels = nothing, duplicates::Symbol = :sum,
)
    rows, samples = _community_record_indices(_community_table_column(table, sample), sample_labels, "sample")
    cols, species_names = _community_record_indices(_community_table_column(table, species), species_labels, "species")
    values = value === nothing ? ones(Int, nrow(table)) : _community_table_column(table, value)
    matrix = _community_triplet_matrix(rows, cols, values, length(samples), length(species_names), duplicates)
    return CommunityData(matrix, samples, species_names, Val(:validated))
end

"""
    align_community(data::CommunityData, tree::CompactTree; missing_tips=:zero)

Return an independent community snapshot with columns in `tree.tip_labels`
order. Unrepresented tree species receive zero columns (`missing_tips=:zero`);
use `missing_tips=:error` to require every tree species in the input. Any input
species absent from the tree causes an error, even if its column is all zero.
The full tree and its root are preserved. Site order, labels, abundances, and
dense/CSC storage are retained. Neither argument is mutated.

Use this once before multiple tree-based calls, or use `prepare_community(data,
tree)` to additionally cache branch sets and occupancy. Direct metric calls
also accept `CommunityData` and align automatically with zero-filled tips.
"""
function align_community(data::CommunityData, tree::CompactTree; missing_tips::Symbol = :zero)
    missing_tips in (:zero, :error) || throw(ArgumentError("missing_tips must be :zero or :error"))
    _validate_community_data(data)
    tree_species = _community_data_labels(tree.tip_labels, tree.ntips, "species")
    tree_set = Set(tree_species)
    unknown = filter(label -> !(label in tree_set), data.species_labels)
    isempty(unknown) || throw(ArgumentError("community species absent from tree: $(join(first(unknown, 5), ", "))"))
    columns = Dict(label => i for (i, label) in enumerate(data.species_labels))
    sources = [get(columns, label, 0) for label in tree_species]
    if missing_tips === :error && any(iszero, sources)
        absent = tree_species[sources .== 0]
        throw(ArgumentError("tree species absent from community: $(join(first(absent, 5), ", "))"))
    end
    matrix = _community_reindex_columns(data.matrix, sources)
    return CommunityData(matrix, copy(data.sample_labels), tree_species, Val(:validated))
end

function _validate_community_data(data::CommunityData)
    _validate_community_values(data.matrix)
    _community_data_labels(data.sample_labels, size(data, 1), "sample")
    _community_data_labels(data.species_labels, size(data, 2), "species")
    return nothing
end

function _community_reindex_columns(matrix::AbstractMatrix{T}, sources) where {T}
    result = zeros(T, size(matrix, 1), length(sources))
    for (target, source) in enumerate(sources)
        source == 0 || copyto!(view(result, :, target), view(matrix, :, source))
    end
    return result
end

function _community_reindex_columns(matrix::SparseMatrixCSC{T}, sources) where {T}
    pointers = ones(Int, length(sources) + 1)
    for (target, source) in enumerate(sources)
        pointers[target + 1] = pointers[target] + (source == 0 ? 0 : length(nzrange(matrix, source)))
    end
    indices = Vector{Int}(undef, last(pointers) - 1)
    values = Vector{T}(undef, length(indices))
    for (target, source) in enumerate(sources)
        source == 0 && continue
        entries = nzrange(matrix, source)
        copyto!(indices, pointers[target], rowvals(matrix), first(entries), length(entries))
        copyto!(values, pointers[target], nonzeros(matrix), first(entries), length(entries))
    end
    return SparseMatrixCSC(size(matrix, 1), length(sources), pointers, indices, values)
end

@inline function _tools_sample_labels(n::Integer, labels::Union{Nothing, AbstractVector{<:AbstractString}})
    labels === nothing && return ["sample_$i" for i in 1:n]
    length(labels) == n || throw(ArgumentError("sample_labels length must match number of sample rows"))
    return String.(labels)
end

function _tools_tip_permutation(tree::CompactTree, tip_labels::Union{Nothing, AbstractVector{<:AbstractString}})
    tip_labels === nothing && return collect(1:tree.ntips)
    length(tip_labels) == tree.ntips || throw(ArgumentError("tip_labels length must match number of sample columns"))
    labels = String.(tip_labels)
    length(unique(labels)) == length(labels) || throw(ArgumentError("tip_labels must not contain duplicates"))
    tip_position = zeros(Int32, tree.nnodes)
    for (position, node) in enumerate(tree.tip_ids)
        tip_position[node] = Int32(position)
    end
    perm = Vector{Int}(undef, tree.ntips)
    for (j, label) in enumerate(labels)
        node = get(tree.tipname_to_id, label, Int32(0))
        node == 0 && throw(ArgumentError("sample column label $label is not present in tree tips"))
        idx = tip_position[node]
        idx == 0 && throw(ArgumentError("tree tip index not found for label $label"))
        perm[j] = idx
    end
    return perm
end

function _tools_tip_permutation(tree::CompactTree, samp::AbstractDataFrame, tip_labels::Nothing)
    labels = String.(names(samp))
    length(labels) == tree.ntips || throw(ArgumentError("dataframe must have one column per tree tip"))
    return _tools_tip_permutation(tree, labels)
end

function _tools_tip_permutation(tree::CompactTree, samp::AbstractDataFrame, tip_labels::AbstractVector{<:AbstractString})
    length(tip_labels) == ncol(samp) || throw(ArgumentError("tip_labels length must match number of sample columns"))
    return _tools_tip_permutation(tree, tip_labels)
end

function _tools_abundance_only(samp::AbstractMatrix{<:Real}, perm::Vector{Int})
    _validate_community_values(samp)
    abund = Matrix{Float64}(undef, size(samp, 1), length(perm))
    @inbounds for i in axes(samp, 1), j in axes(samp, 2)
        abund[i, perm[j]] = Float64(samp[i, j])
    end
    return abund
end

function _tools_presence_only(samp::AbstractMatrix{<:Real}, perm::Vector{Int})
    _validate_community_values(samp)
    presence = Matrix{Bool}(undef, size(samp, 1), length(perm))
    @inbounds for i in axes(samp, 1), j in axes(samp, 2)
        presence[i, perm[j]] = samp[i, j] > 0
    end
    return presence
end

# Legacy consumers need dense workspaces, but CSC inputs only need stored-entry
# traversal. SES bypasses these conversions and uses sparse presence rows.
function _tools_abundance_only(samp::SparseMatrixCSC{<:Real}, perm::Vector{Int})
    _validate_community_values(samp)
    abundance = zeros(Float64, size(samp, 1), length(perm))
    for col in axes(samp, 2), entry in nzrange(samp, col)
        abundance[rowvals(samp)[entry], perm[col]] = Float64(nonzeros(samp)[entry])
    end
    return abundance
end

function _tools_presence_only(samp::SparseMatrixCSC{<:Real}, perm::Vector{Int})
    _validate_community_values(samp)
    presence = fill(false, size(samp, 1), length(perm))
    for col in axes(samp, 2), entry in nzrange(samp, col)
        presence[rowvals(samp)[entry], perm[col]] = nonzeros(samp)[entry] > 0
    end
    return presence
end

function _tools_presence_matrix(samp::AbstractMatrix{<:Real}, perm::Vector{Int})
    return _tools_presence_only(samp, perm)
end

function _community_abundance(samp::AbstractMatrix{<:Real}, tree::CompactTree, tip_labels)
    size(samp, 2) == tree.ntips || throw(ArgumentError("community matrix must have one column per tree tip"))
    perm = _tools_tip_permutation(tree, tip_labels)
    return _tools_abundance_only(samp, perm)
end

function _community_presence(samp::AbstractMatrix{<:Real}, tree::CompactTree, tip_labels)
    size(samp, 2) == tree.ntips || throw(ArgumentError("community matrix must have one column per tree tip"))
    perm = _tools_tip_permutation(tree, tip_labels)
    return _tools_presence_only(samp, perm)
end

# The presence-based kernels below visit stored entries of CSC input directly.
# Dense conversion helpers above serve older diversity/distance consumers;
# SES now uses sparse rows and does not call those conversion helpers.
function _validate_community_values(community::AbstractMatrix{<:Real})
    Base.require_one_based_indexing(community)
    values = community isa SparseMatrixCSC ? nonzeros(community) : community
    all(x -> isfinite(x) && x >= 0, values) ||
        throw(ArgumentError("community values must be finite and non-negative"))
    return nothing
end

function _validate_community_branch_lengths(tree::CompactTree)
    all(x -> isfinite(x) && x >= 0, tree.edge_length) ||
        throw(ArgumentError("community metrics require finite, non-negative branch lengths"))
    isfinite(sum(tree.edge_length)) ||
        throw(ArgumentError("total branch length must be finite"))
    return nothing
end

function _community_presence_rows(community::AbstractMatrix{<:Real}, perm::AbstractVector{<:Integer})
    _validate_community_values(community)
    rows = [Int32[] for _ in axes(community, 1)]
    for col in axes(community, 2), site in axes(community, 1)
        community[site, col] > 0 && push!(rows[site], Int32(perm[col]))
    end
    return rows
end

function _community_presence_rows(community::SparseMatrixCSC{<:Real}, perm::AbstractVector{<:Integer})
    _validate_community_values(community)
    rows = [Int32[] for _ in axes(community, 1)]
    row_indices, values = rowvals(community), nonzeros(community)
    for col in axes(community, 2), entry in nzrange(community, col)
        values[entry] > 0 && push!(rows[row_indices[entry]], Int32(perm[col]))
    end
    return rows
end

function _community_tip_rows(community::AbstractMatrix{<:Real}, tree::CompactTree, tip_labels)
    size(community, 2) == tree.ntips ||
        throw(ArgumentError("community matrix must have one column per tree tip"))
    perm = _tools_tip_permutation(tree, tip_labels)
    return _community_tip_rows(community, perm)
end

function _community_tip_rows(community::AbstractMatrix{<:Real}, perm::AbstractVector{<:Integer})
    rows = _community_presence_rows(community, perm)
    # Match tree-tip order even when input columns were permuted.
    if !issorted(perm)
        foreach(sort!, rows)
    end
    return rows
end
