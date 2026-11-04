"""
    TreeSampleResult

Iterable sample from [`tree_sampler`](@ref). `trees`, `sample_indices`,
`source_paths`, `source_indices` and `source_lines` follow random draw order.
Sample indices address the eligible population concatenated across files;
source indices are original one-based record numbers within each file.
`population_size` counts eligible records; `replace`, `burnin` and `thin`
record the sampling options. Trees and metadata arrays are independently owned.
"""
struct TreeSampleResult
    trees::Vector{CompactTree}
    sample_indices::Vector{Int}
    source_paths::Vector{String}
    source_indices::Vector{Int}
    source_lines::Vector{Int}
    population_size::Int
    replace::Bool
    burnin::Int
    thin::Int
end

Base.length(result::TreeSampleResult) = length(result.trees)
Base.eltype(::Type{TreeSampleResult}) = CompactTree
Base.iterate(result::TreeSampleResult, state=1) =
    state > length(result) ? nothing : (result.trees[state], state + 1)

function _tree_sample_indices(rng, population, n, replace)
    replace && return [rand(rng, 1:population) for _ in 1:n]
    n <= population || throw(ArgumentError("n exceeds the eligible population for sampling without replacement"))
    # Partial Fisher-Yates over a virtual index vector: O(n) storage even when
    # sampling a few trees from a very large posterior.
    moved = Dict{Int, Int}()
    indices = Vector{Int}(undef, n)
    for i in 1:n
        remaining = population - i + 1
        selected = rand(rng, 1:remaining)
        indices[i] = get(moved, selected, selected)
        moved[selected] = get(moved, remaining, remaining)
        delete!(moved, remaining)
    end
    return indices
end

function _tree_file_stamps(paths)
    return map(paths) do path
        info = stat(path)
        (; info.size, info.mtime)
    end
end

"""
    tree_sampler(paths; n=100, replace=true, burnin=0, thin=1, rng=Random.default_rng())

Sample posterior trees uniformly from one or several Newick files, returning
an iterable `TreeSampleResult`. `paths`, record syntax, and per-file burn-in /
thinning follow [`read_trees`](@ref). Each eligible tree has equal probability,
including when files contain different numbers of trees. The default draws
with replacement, matching phyloregion's `tree_sampler`; `replace=false` draws
distinct records. `n` must be a positive integer. Empty eligible populations
and requesting more distinct trees than available are errors.

Use an explicit `rng` for reproducible draw indices; Julia and R do not share
random sequences. Returned trees follow draw order. Repeated draws are parsed
once then copied, so modifying one does not change another. Provenance retains
the original file, record number, line, and eligible-population index.

Two sequential scans avoid storing all records: the first counts eligible
records and checks framing; the second parses only sampled records. Structural
or branch-token errors in unsampled/burned/thinned records are not parsed.
Unterminated records, quotes and comments still fail during scanning. All
files must remain unchanged between scans; size/time and population changes
are checked. File handles close on completion or failure.

Memory is O(n + current record + returned trees), plus source paths; sampling
without replacement does not allocate a vector the size of the population.
Pass the result to a generator such as
`mean_dist(phylobeta(data, tree) for tree in sample)` to accumulate distances
without retaining every distance matrix. Use `read_trees` to process all
eligible trees with only one parsed tree at a time.
"""
function tree_sampler(paths; n=100, replace::Bool=true, burnin=0, thin=1,
    rng::AbstractRNG=Random.default_rng(),
)
    count = _tree_stream_integer(n, "n", 1)
    stream = read_trees(paths; burnin, thin)
    stamps = _tree_file_stamps(stream.paths)
    population = 0
    try
        while _next_tree_record!(stream) !== nothing
            population += 1
        end
    finally
        close(stream)
    end
    population > 0 || throw(ArgumentError("no eligible trees remain after burnin and thinning"))
    indices = _tree_sample_indices(rng, population, count, replace)
    stamps == _tree_file_stamps(stream.paths) || throw(ArgumentError("tree files changed during sampling"))
    order = sortperm(indices)
    trees = Vector{CompactTree}(undef, count)
    source_paths, source_indices, source_lines = Vector{String}(undef, count), zeros(Int, count), zeros(Int, count)
    position, seen = 1, 0
    read_trees(stream.paths; burnin=stream.burnin, thin=stream.thin) do source
        while (record = _next_tree_record!(source)) !== nothing
            seen += 1
            position <= count && indices[order[position]] == seen || continue
            tree = _parse_tree_record(record)
            first_draw = true
            while position <= count && indices[order[position]] == seen
                draw = order[position]
                trees[draw] = first_draw ? tree : deepcopy(tree)
                source_paths[draw], source_indices[draw], source_lines[draw] = record.path, record.index, record.line
                first_draw = false
                position += 1
            end
        end
    end
    seen == population && position == count + 1 && stamps == _tree_file_stamps(stream.paths) ||
        throw(ArgumentError("tree files changed between sampling scans"))
    return TreeSampleResult(trees, indices, source_paths, source_indices, source_lines,
        population, replace, stream.burnin, stream.thin)
end
