"""
    PhylobetaSESResult

Standardized effect size of one phylogenetic beta component, returned by
[`phylobeta_ses`](@ref). `observed`, `null_mean`, `null_sd`, `z_score`, `p_value`
and `valid_replicates` are symmetric matrices, or vectors in requested `pairs`
order. `richness` and `sample_labels` describe the complete input site order.
Arrays are independently owned; undefined values are NaN, not imputed zeros.
"""
struct PhylobetaSESResult{N}
    sample_labels::Vector{String}
    richness::Vector{Int}
    observed::Array{Float64, N}
    null_mean::Array{Float64, N}
    null_sd::Array{Float64, N}
    z_score::Array{Float64, N}
    p_value::Array{Float64, N}
    valid_replicates::Array{Int, N}
    pairs::Union{Nothing, Vector{Tuple{Int, Int}}}
    index_family::Symbol
    component::Symbol
    null_model::Symbol
    runs::Int
    iterations::Int
    max_attempts_per_swap::Int
end

# Reuse branch vectors and one ancestor marker across sites and replicates.
# For selected pairs only their endpoints need branch sets; randomization still
# uses the complete community so global frequency/swap constraints are retained.
mutable struct _BetaSESWorkspace
    branches::Vector{Vector{Int32}}
    sites::Vector{Int}
    edge_of_child::Vector{Int32}
    seen::Vector{Int}
    generation::Int
end

function _BetaSESWorkspace(tree, n, pairs)
    used = trues(n)
    if pairs !== nothing
        fill!(used, false)
        for (i, j) in pairs
            used[i] = used[j] = true
        end
    end
    return _BetaSESWorkspace([Int32[] for _ in 1:n], findall(used),
        _community_edge_index(tree), zeros(Int, tree.nnodes), 0)
end

function _beta_ses_branches!(workspace, tree, rows)
    for site in workspace.sites
        workspace.generation == typemax(Int) && (fill!(workspace.seen, 0); workspace.generation = 0)
        workspace.generation += 1
        edges = workspace.branches[site]
        empty!(edges)
        for tip in rows[site]
            node = tree.tip_ids[tip]
            while node != tree.root && workspace.seen[node] != workspace.generation
                workspace.seen[node] = workspace.generation
                push!(edges, workspace.edge_of_child[node])
                node = tree.parent_of_node[node]
            end
        end
        sort!(edges)
    end
    return workspace.branches
end

function _beta_ses_visit(f, workspace, tree, family, component, pairs)
    branches = workspace.branches
    if pairs === nothing
        n = length(branches)
        for j in 1:n, i in j:n
            abc = _community_branch_components(branches[i], branches[j], tree.edge_length)
            f(i + (j - 1) * n, _community_beta_components(abc..., family)[component])
        end
    else
        for (index, (i, j)) in enumerate(pairs)
            # Use the same operand order as the full lower triangle. Otherwise
            # a reversed pair can round differently and change exact rank ties.
            i, j = i < j ? (j, i) : (i, j)
            abc = _community_branch_components(branches[i], branches[j], tree.edge_length)
            f(index, _community_beta_components(abc..., family)[component])
        end
    end
    return nothing
end

function _beta_ses_parameter(value, name, minimum)
    value isa Integer && !(value isa Bool) && minimum <= value <= typemax(Int) ||
        throw(ArgumentError("$name must be an integer between $minimum and typemax(Int)"))
    return Int(value)
end

function _beta_ses_model(model)
    model === :tipshuffle && return :taxa_labels
    model in (:rowwise, :colwise) && throw(ArgumentError(
        "R rowwise/colwise permutations are undone by label alignment; use :frequency or :richness for actual occurrence randomization"))
    return _tools_normalize_phylo_null(model)
end

"""
    phylobeta_ses(community, tree; index_family=:sorensen, component=:total,
        null_model=:taxa_labels, runs=999, pairs=nothing, iterations=1000,
        max_attempts_per_swap=10_000, tip_labels=nothing, sample_labels=nothing,
        rng=Random.default_rng())

Standardize one rooted phylogenetic beta component against a sparse null model.
Select :total, :turnover or :nestedness and :sorensen or :jaccard. Returns
`PhylobetaSESResult` with observed values, finite-draw mean and sample SD,
z=(observed-mean)/SD, lower-tail mid-rank p and finite replicate counts. To
compare components under the same draws, use fresh RNGs with the same seed.

Models share the `ses_pd` sparse engine: :taxa_labels (alias :tipshuffle) uses
one common permutation of the complete tree-tip pool, preserving all taxonomic
pairwise overlap counts. :richness and :sample_pool independently select each
site's observed richness from that complete pool; :phylogeny_pool additionally
applies a common label permutation. :frequency preserves each species' site
occupancy and randomizes its sites. :trialswap/:independentswap perform bounded
checkerboard trials preserving both margins, starting from the observed matrix
for every replicate. These finite swap runs need not mix to a uniform null.
R's :rowwise/:colwise reorder labeled rows/columns and are undone by alignment;
they are rejected rather than treated as new ecological null models.

Dense and CSC input values must be finite/non-negative; positive abundance
denotes presence. Matrix inputs require one column per tree tip, including
zero columns. `CommunityData` aligns to the complete tree and fills absent
species with zeros. All tree tips remain eligible, including unobserved ones.
No pruning, silent label intersection or name rewriting is performed. Site
labels are unique nonempty strings/integers. Inputs and the tree are unchanged.

Full output is symmetric in input site order. `pairs` returns only requested
vectors, retaining repeated, reversed and self pairs. Randomization always
uses all input sites/species, including sites absent from the requested pairs.
Only branch sets of requested endpoints are constructed. Working memory scales
with presence/branch sets, tree size and output size, never sites*species or
output pairs*runs. Online Welford statistics retain no replicate matrices.

Each pair uses only finite null values. With zero valid draws, null mean/SD/z/p
are NaN; with one, SD/z are NaN. Constant nulls retain SD=0 and yield NaN z.
Nonfinite observed values also yield NaN z/p, even when their null is defined.
Finite observed p=(1+number(null<observed)+0.5number(null==observed))/(valid+1).
This is a lower-tail rank score with exact ties, not a two-sided p-value;
its denominator uses valid draws, unlike the older site-SES runs+1 convention.
Empty/zero-PD sites follow `phylobeta` component rules, including NaN diagonals.
SES z scores can be negative/undefined and are not distances for regionalization.
`rng` controls all draws; fixed seed, input order and null parameters reproduce
the sequence. Requested pairs and component do not alter RNG consumption.
"""
function phylobeta_ses(community::AbstractMatrix{<:Real}, tree::CompactTree;
    index_family::Symbol = :sorensen, component::Symbol = :total,
    null_model::Symbol = :taxa_labels, runs = 999, pairs = nothing,
    iterations = 1000, max_attempts_per_swap = 10_000,
    tip_labels::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    sample_labels = nothing, rng::AbstractRNG = Random.default_rng(),
)
    _validate_community_beta_family(index_family)
    component_index = findfirst(==(component), (:total, :turnover, :nestedness))
    component_index === nothing && throw(ArgumentError("component must be :total, :turnover or :nestedness"))
    count_runs = _beta_ses_parameter(runs, "runs", 1)
    count_iterations = _beta_ses_parameter(iterations, "iterations", 0)
    attempts = _beta_ses_parameter(max_attempts_per_swap, "max_attempts_per_swap", 1)
    model = _beta_ses_model(null_model)
    _validate_community_branch_lengths(tree)
    n = size(community, 1)
    labels = _validated_labels(sample_labels, n, "sample")
    selected = _community_pairs(pairs, n)
    rows = _community_tip_rows(community, tree, tip_labels)
    null = _CommunityNullWorkspace(rows, tree.ntips)
    workspace = _BetaSESWorkspace(tree, n, selected)
    shape = selected === nothing ? (n, n) : (length(selected),)
    observed = fill(NaN, shape)
    stats = _SESOnlineStats(shape)
    _beta_ses_branches!(workspace, tree, rows)
    _beta_ses_visit(workspace, tree, index_family, component_index, selected) do i, value
        observed[i] = value
    end
    for _ in 1:count_runs
        random = _randomize_community!(null, rows, model, count_iterations, attempts, rng)
        _beta_ses_branches!(workspace, tree, random)
        _beta_ses_visit(workspace, tree, index_family, component_index, selected) do i, value
            _ses_online_update!(stats, i, observed[i], value)
        end
    end
    indices = selected === nothing ? (i + (j - 1) * n for j in 1:n for i in j:n) : eachindex(observed)
    (; means, sds, z, ranks, counts) = _ses_online_finish!(stats, observed, indices)
    p = fill(NaN, shape)
    for i in indices
        if isfinite(observed[i]) && counts[i] > 0
            p[i] = ranks[i] / (counts[i] + 1.0)
        end
    end
    if selected === nothing
        for values in (observed, means, sds, z, p, counts), j in 1:n, i in (j + 1):n
            values[j, i] = values[i, j]
        end
    end
    return PhylobetaSESResult(labels, copy(null.richness), observed, means, sds, z, p,
        counts, selected, index_family, component, model, count_runs, count_iterations, attempts)
end

function phylobeta_ses(data::CommunityData, tree::CompactTree; kwargs...)
    _community_data_keywords(kwargs)
    aligned = _community_tree_data(data, tree)
    return phylobeta_ses(aligned.matrix, tree; sample_labels = aligned.sample_labels, kwargs...)
end
