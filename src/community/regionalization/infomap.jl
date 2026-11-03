"""
    InfomapResult

Network regionalization returned by [`infomap`](@ref). `sample_labels` and
`membership` retain input row order. Regions `1:k` are numbered by their first
site, and `cluster_sizes` counts sites in that order. `codelength` is igraph's
reported code length in bits for its complete site–species network partition,
including species nodes; it is not a score for the site labels alone.

The result records `weights`, `trials`, `seed`, `regularized`,
`regularization_strength`, and the native `igraph_version`. Arrays are owned by
the result. There is no hierarchy or distance matrix in this result.
"""
struct InfomapResult
    sample_labels::Vector{String}
    membership::Vector{Int}
    cluster_sizes::Vector{Int}
    k::Int
    codelength::Float64
    weights::Symbol
    trials::Int
    seed::UInt64
    regularized::Bool
    regularization_strength::Float64
    igraph_version::String
end

function _infomap_settings(weights, trials, seed, regularized, strength)
    weights in (:binary, :abundance) ||
        throw(ArgumentError("weights must be :binary or :abundance"))
    trials isa Integer && !(trials isa Bool) && 1 <= trials <= typemax(Int) ||
        throw(ArgumentError("trials must be a positive Int-sized integer"))
    seed isa Integer && !(seed isa Bool) && 0 <= seed <= typemax(UInt64) ||
        throw(ArgumentError("seed must be an integer in 0:typemax(UInt64)"))
    regularized isa Bool || throw(ArgumentError("regularized must be a Bool"))
    strength isa Real && !(strength isa Bool) && isfinite(strength) && strength > 0 ||
        throw(ArgumentError("regularization_strength must be finite and positive"))
    converted = Float64(strength)
    isfinite(converted) && converted > 0 ||
        throw(ArgumentError("regularization_strength must be representable as a positive finite Float64"))
    return Int(trials), UInt64(seed), converted
end

# CSC visits only stored entries (including explicit zeros); dense inputs are
# scanned once. Unobserved species are omitted without changing site order.
# Endpoint IDs are zero-based for igraph, with disjoint site/species ranges.
function _infomap_network(community::AbstractMatrix{<:Real}, weights)
    n, p = size(community)
    n > 0 && p > 0 || throw(ArgumentError("Infomap requires sites and species"))
    _validate_community_values(community)
    species_nodes = zeros(Int, p)
    active_species = Int[]
    richness = zeros(Int, n)
    edges = Int[]
    edge_weights = weights === :abundance ? Float64[] : nothing
    total_weight = 0.0
    _community_each_positive(community) do row, col, value
        node = species_nodes[col]
        if node == 0
            push!(active_species, col)
            node = Base.checked_add(n, length(active_species))
            species_nodes[col] = node
        end
        push!(edges, row - 1, node - 1)
        richness[row] += 1
        if edge_weights !== nothing
            converted = Float64(value)
            isfinite(converted) && converted > 0 ||
                throw(ArgumentError("abundances must be representable as positive finite Float64 weights"))
            total_weight += converted
            isfinite(2 * total_weight) ||
                throw(ArgumentError("undirected weight sums exceed finite Float64 range; rescale abundances"))
            push!(edge_weights, converted)
        end
    end
    all(>(0), richness) ||
        throw(ArgumentError("Infomap does not assign empty sites; select nonempty sites and their labels explicitly"))
    return (; edges, edge_weights, nvertices = Base.checked_add(n, length(active_species)),
        active_species, richness)
end

function _infomap_site_membership(nodes, n)
    membership = Vector{Int}(undef, n)
    regions = Dict{Int, Int}()
    sizes = Int[]
    for site in 1:n
        node = nodes[site]
        region = get(regions, node, 0)
        if region == 0
            push!(sizes, 0)
            region = length(sizes)
            regions[node] = region
        end
        membership[site] = region
        sizes[region] += 1
    end
    return membership, sizes
end

"""
    infomap(community; weights=:binary, trials=10, seed=20260713,
            regularized=true, regularization_strength=1.0,
            sample_labels=nothing, species_labels=nothing)
    infomap(data::CommunityData; ...)

Find site regions with igraph Infomap on an undirected site–species bipartite
network. Rows are sites, columns are species; values must be finite and
nonnegative. `weights=:binary` makes each positive entry one edge, independent
of its magnitude. `:abundance` uses positive finite Float64 abundances as edge
weights, with finite undirected strength sums. CSC input is processed through
its stored entries without a dense community, adjacency or beta matrix.

Regularization is enabled by default. `regularization_strength` is finite and
positive; it is ignored by the engine when `regularized=false`. Regularization
can change both memberships and the number of regions, even on disconnected
networks. It does not enforce bipartite restrictions on its prior or geographic
contiguity. Input zeroes mean no observed edge, not missing survey effort.

Empty sites are rejected. Unobserved species columns are omitted before
fitting; explicit stored zeroes do not create edges. Labels are validated
separately for each node type, so a site and a species may share a name.
`CommunityData` owns its labels; set them when constructing it. For lineage or
functional regionalization, pass the aggregated `CommunityData` returned by
`collapse_range` or `collapse_traits`.

The engine searches `trials` times and reports its best partition. A private
native PCG32 RNG uses `seed` (0 through typemax(UInt64)); Julia's RNG and the
previous native RNG state are preserved. EvoKits serializes its native calls.
Identical data order, settings and engine version are reproducible, but
different engines or input permutations need not select the same local
optimum. Neither a trial count nor regularization proves ecological validity.
Inputs are not mutated. Returns an [`InfomapResult`](@ref).
"""
function infomap(community::AbstractMatrix{<:Real}; weights=:binary, trials=10,
    seed=20260713, regularized=true, regularization_strength=1.0,
    sample_labels=nothing, species_labels=nothing)
    attempts, native_seed, strength = _infomap_settings(weights, trials, seed,
        regularized, regularization_strength)
    labels = _validated_labels(sample_labels, size(community, 1), "sample")
    _validated_labels(species_labels, size(community, 2), "species")
    network = _infomap_network(community, weights)
    nodes, codelength, version = _infomap_backend(network.edges, network.edge_weights,
        network.nvertices, attempts, native_seed, regularized, strength)
    membership, sizes = _infomap_site_membership(nodes, size(community, 1))
    return InfomapResult(labels, membership, sizes, length(sizes), codelength,
        weights, attempts, native_seed, regularized, strength, version)
end

function infomap(data::CommunityData; kwargs...)
    _community_aggregation_keywords(kwargs)
    return infomap(data.matrix; sample_labels=data.sample_labels,
        species_labels=data.species_labels, kwargs...)
end

infomap(community::AbstractMatrix; kwargs...) =
    throw(ArgumentError("community matrix must contain real numbers without missing values"))
