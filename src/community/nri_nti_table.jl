"""Flat, sorted tip-position communities read from a contiguous-block CSV table."""
struct _NriTableCommunityData
    grid_ids::Vector{String}
    tip_positions::Vector{Int32}
    offsets::Vector{Int64}
end

_nri_table_view(data::_NriTableCommunityData, index::Int) =
    @view data.tip_positions[data.offsets[index]:(data.offsets[index + 1] - 1)]

function _nri_tip_position_index(tree::CompactTree)
    positions = zeros(Int32, tree.nnodes)
    for (position, node) in enumerate(tree.tip_ids)
        positions[Int(node)] = Int32(position)
    end
    return positions
end

function _finish_nri_table_grid!(grid_ids::Vector{String}, tip_positions::Vector{Int32},
                               offsets::Vector{Int64}, completed::Set{String},
                               grid_id::String, current::Vector{Int32})
    sort!(current)
    unique!(current)
    append!(tip_positions, current)
    push!(offsets, Int64(length(tip_positions) + 1))
    push!(grid_ids, grid_id)
    push!(completed, grid_id)
    empty!(current)
    return nothing
end

function _read_nri_table_communities(tree::CompactTree, path::AbstractString,
                                     dataset::AbstractString, delim::Union{Nothing,Char})
    isfile(path) || throw(ArgumentError("$dataset table does not exist: $path"))
    rows = CSV.Rows(path; reusebuffer = true, delim = delim)

    tip_position_of_node = _nri_tip_position_index(tree)
    grid_ids = String[]
    tip_positions = Int32[]
    offsets = Int64[1]
    completed = Set{String}()
    current = Int32[]
    current_grid = ""
    have_grid = false
    row_number = 0

    for row in rows
        row_number += 1
        hasproperty(row, :grid_id) ||
            throw(ArgumentError("$dataset table must contain a grid_id column"))
        hasproperty(row, :species) ||
            throw(ArgumentError("$dataset table must contain a species column"))
        ismissing(row.grid_id) && throw(ArgumentError("$dataset table has a missing grid_id at row $row_number"))
        ismissing(row.species) && throw(ArgumentError("$dataset table has a missing species at row $row_number"))
        grid_id = String(row.grid_id)
        species = String(row.species)
        if !have_grid || grid_id != current_grid
            if have_grid
                _finish_nri_table_grid!(grid_ids, tip_positions, offsets, completed,
                                        current_grid, current)
            end
            grid_id in completed &&
                throw(ArgumentError("$dataset table grid '$grid_id' occurs in non-contiguous blocks"))
            current_grid = grid_id
            have_grid = true
        end
        node = get(tree.tipname_to_id, species, Int32(0))
        node == 0 && throw(ArgumentError("$dataset table species '$species' is not a tree tip"))
        push!(current, tip_position_of_node[Int(node)])
    end
    have_grid && _finish_nri_table_grid!(grid_ids, tip_positions, offsets, completed,
                                         current_grid, current)
    return _NriTableCommunityData(grid_ids, tip_positions, offsets)
end

function _nri_table_subset(pool, observed)
    pool_index = firstindex(pool)
    pool_last = lastindex(pool)
    for position in observed
        while pool_index <= pool_last && pool[pool_index] < position
            pool_index += 1
        end
        pool_index <= pool_last && pool[pool_index] == position || return false
    end
    return true
end

"""
    nri_nti_table(tree, pool_table, observed_table; nreps = 999, seed = 20260713,
                  threaded = true, min_observed_richness = 2, delim = nothing)

Run local-pool NRI/NTI for every eligible grid in two delimited text tables.
Each table must have `grid_id,species` columns and each grid must occupy one
contiguous block. `delim = nothing` lets CSV.jl detect the delimiter; use, for
example, `delim = '\t'` for TSV. The return value has `grid_ids` and same-order
`results`. When `threaded = true`, grids are dynamically distributed across
Julia threads; per-grid seeds are stable regardless of execution order.
"""
function nri_nti_table(
    tree::CompactTree,
    pool_table::AbstractString,
    observed_table::AbstractString;
    nreps::Integer = 999,
    seed::Integer = 20260713,
    threaded::Bool = true,
    min_observed_richness::Integer = 2,
    delim::Union{Nothing,Char} = nothing,
)
    nreps >= 2 || throw(ArgumentError("nreps must be at least 2 to calculate sample standard deviations"))
    min_observed_richness >= 2 || throw(ArgumentError("min_observed_richness must be at least 2"))
    all(isfinite, tree.edge_length) && all(tree.edge_length .>= 0.0) ||
        throw(ArgumentError("tree contains invalid branch lengths"))

    pool_data = _read_nri_table_communities(tree, pool_table, "pool", delim)
    observed_data = _read_nri_table_communities(tree, observed_table, "observed", delim)
    pool_index = Dict(grid_id => i for (i, grid_id) in enumerate(pool_data.grid_ids))
    grid_ids = String[]
    pool_indices = Int[]
    observed_indices = Int[]
    for (observed_index, grid_id) in enumerate(observed_data.grid_ids)
        pool_position = get(pool_index, grid_id, 0)
        pool_position == 0 && throw(ArgumentError("observed grid '$grid_id' has no local pool"))
        pool = _nri_table_view(pool_data, pool_position)
        observed = _nri_table_view(observed_data, observed_index)
        _nri_table_subset(pool, observed) ||
            throw(ArgumentError("observed species are not a subset of the local pool for grid '$grid_id'"))
        length(observed) >= min_observed_richness || continue
        push!(grid_ids, grid_id)
        push!(pool_indices, pool_position)
        push!(observed_indices, observed_index)
    end

    query_index = _build_tree_query_index(tree)
    workspaces = [_NriMetricsWorkspace(tree) for _ in 1:Threads.maxthreadid()]
    results = Vector{NriNtiResult}(undef, length(grid_ids))
    function run_grid!(index::Int)
        grid_id = grid_ids[index]
        pool = _nri_table_view(pool_data, pool_indices[index])
        observed = _nri_table_view(observed_data, observed_indices[index])
        effective_seed = stable_grid_seed(seed, grid_id)
        results[index] = _nri_nti_positions(tree, pool, observed, query_index,
                                             workspaces[Threads.threadid()];
                                             nreps, random=StableRNG(effective_seed), effective_seed)
        return nothing
    end

    if threaded && Threads.nthreads() > 1
        Threads.@threads :dynamic for index in eachindex(grid_ids)
            run_grid!(index)
        end
    else
        for index in eachindex(grid_ids)
            run_grid!(index)
        end
    end
    return (grid_ids=grid_ids, results=results)
end
