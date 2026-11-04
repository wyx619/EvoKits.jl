"""
    community_longtable(data::CommunityData)

Export a labeled community as a named tuple with `table`, `sample_labels`, and
`species_labels`. The `DataFrame` has `sample`, `species`, and `value` columns,
with one row per positive site/species combination, ordered by species column
and then site row. Values retain their numeric type and abundance, including
fractional weights; abundances are not expanded into repeated occurrences.
Zero values, including explicit CSC zeros, are omitted.

The complete label lists preserve the original order, empty sites and all-zero
species columns. Keep them with the table when saving or reconstructing data:

```julia
records = community_longtable(data)
restored = CommunityData(records.table;
    sample_labels=records.sample_labels, species_labels=records.species_labels)
```

Reconstruction uses CSC storage and the usual long-table value rules (Bool
becomes integer counts). The table alone cannot preserve empty label domains
or the original site order. Output arrays are independent of the input. CSC
input is traversed directly, with storage proportional to the positive records
and complete label lists, without building a dense matrix or occurrence table.
"""
function community_longtable(data::CommunityData)
    _validate_community_data(data)
    matrix = data.matrix
    values = matrix isa SparseMatrixCSC ? nonzeros(matrix) : matrix
    n = count(>(0), values)
    samples = Vector{String}(undef, n)
    species = Vector{String}(undef, n)
    abundances = Vector{eltype(matrix)}(undef, n)
    cursor = Ref(0)
    _community_each_positive(matrix) do row, col, value
        index = (cursor[] += 1)
        samples[index] = data.sample_labels[row]
        species[index] = data.species_labels[col]
        abundances[index] = value
    end
    table = DataFrame(sample = samples, species = species, value = abundances;
        copycols = false)
    return (; table, sample_labels = copy(data.sample_labels),
        species_labels = copy(data.species_labels))
end
