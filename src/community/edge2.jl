"""
    Edge2Result

Result of an EDGE2 calculation. The vectors follow `tree.tip_labels` order.
"""
Base.@kwdef struct Edge2Result
    species::Vector{String} = String[]
    tbl::Vector{Float64} = Float64[]
    pext::Vector{Float64} = Float64[]
    ed::Vector{Float64} = Float64[]
    edge::Vector{Float64} = Float64[]
    pd::Float64 = 0.0
    epd_loss::Float64 = 0.0
end

# Keep the field names used by the original reference module available while
# exposing idiomatic lower-case fields to new EvoKits callers.
function Base.getproperty(result::Edge2Result, name::Symbol)
    name === :TBL && return getfield(result, :tbl)
    name === :ED && return getfield(result, :ed)
    name === :EDGE && return getfield(result, :edge)
    name === :PD && return getfield(result, :pd)
    name === :ePDloss && return getfield(result, :epd_loss)
    return getfield(result, name)
end

struct PextDistribution
    categories::Vector{String}
    pext_values::Vector{Float64}
    function PextDistribution(categories::AbstractVector{<:AbstractString}, pext_values::AbstractVector{<:Real})
        length(categories) == length(pext_values) || throw(ArgumentError("categories and pext_values must have the same length"))
        isempty(categories) && throw(ArgumentError("pext distribution cannot be empty"))
        values = Float64.(pext_values)
        all(isfinite, values) && all((0.0 .< values) .& (values .<= 1.0)) ||
            throw(ArgumentError("pext_values must be finite and in (0, 1]"))
        new(String.(categories), values)
    end
end

function _edge2_validate_pext(tree::CompactTree, species_pext::AbstractVector{<:Real})
    length(species_pext) == tree.ntips || throw(ArgumentError("species_pext must have one value per tree tip"))
    pext = Float64.(species_pext)
    all(isfinite, pext) && all((0.0 .< pext) .& (pext .<= 1.0)) ||
        throw(ArgumentError("species_pext must be finite and in (0, 1]"))
    @inbounds for edge in 1:tree.nedges
        length = tree.edge_length[edge]
        isfinite(length) && length >= 0.0 || throw(ArgumentError("tree contains invalid branch lengths"))
    end
    return pext
end

"""
    edge2_mod(tree::CompactTree, species_pext) -> Edge2Result

Compute Evolutionary Distinct and Globally Endangered (EDGE2) scores. The
extinction probabilities must follow `tree.tip_labels` order.
"""
function edge2_mod(tree::CompactTree, species_pext::AbstractVector{<:Real})
    pext = _edge2_validate_pext(tree, species_pext)
    tree.ntips > 0 || throw(ArgumentError("tree must contain at least one tip"))

    node_pext = ones(Float64, tree.nnodes)
    @inbounds for (tip_position, tip) in enumerate(tree.tip_ids)
        node_pext[Int(tip)] = pext[tip_position]
    end
    @inbounds for node in tree.postorder
        inode = Int(node)
        inode == Int(tree.root) && continue
        node_pext[Int(tree.parent_of_node[inode])] *= node_pext[inode]
    end

    weighted_edge = zeros(Float64, tree.nedges)
    edge_of_child = zeros(Int32, tree.nnodes)
    @inbounds for edge in 1:tree.nedges
        child = Int(tree.child_of_edge[edge])
        weighted_edge[edge] = tree.edge_length[edge] * node_pext[child]
        edge_of_child[child] = Int32(edge)
    end

    ancestral_edge_sum = zeros(Float64, tree.nnodes)
    @inbounds for node in tree.preorder
        inode = Int(node)
        inode == Int(tree.root) && continue
        parent = Int(tree.parent_of_node[inode])
        edge = Int(edge_of_child[inode])
        ancestral_edge_sum[inode] = ancestral_edge_sum[parent] + weighted_edge[edge]
    end

    tbl = Vector{Float64}(undef, tree.ntips)
    edge = Vector{Float64}(undef, tree.ntips)
    @inbounds for (tip_position, tip) in enumerate(tree.tip_ids)
        inode = Int(tip)
        terminal_edge = Int(edge_of_child[inode])
        tbl[tip_position] = tree.edge_length[terminal_edge]
        edge[tip_position] = ancestral_edge_sum[Int(tree.parent_of_node[inode])] + weighted_edge[terminal_edge]
    end

    return Edge2Result(
        species = copy(tree.tip_labels),
        tbl = tbl,
        pext = pext,
        ed = edge ./ pext,
        edge = edge,
        pd = sum(tree.edge_length),
        epd_loss = sum(weighted_edge),
    )
end

function _edge2_category_pool(distribution::PextDistribution, categories::AbstractVector{<:AbstractString})
    mask = falses(length(distribution.categories))
    for category in categories
        mask .|= distribution.categories .== String(category)
    end
    pool = distribution.pext_values[mask]
    isempty(pool) && throw(ArgumentError("pext distribution has no values for categories $(collect(categories))"))
    return pool
end

"""
    calc_pext_distribution(; anchor_pext, n_tips, rng) -> PextDistribution

Generate a sampled extinction-probability distribution using the EDGE2
category calibration procedure.
"""
function calc_pext_distribution(; anchor_pext::AbstractVector{<:Real} = [0.97, 0.485, 0.2425, 0.12125, 0.060625], n_tips::Integer = 10_000, rng::AbstractRNG = Random.default_rng())
    length(anchor_pext) == 5 || throw(ArgumentError("anchor_pext must contain five category anchors"))
    n_tips >= 5 || throw(ArgumentError("n_tips must be at least 5"))
    anchors = Float64.(anchor_pext)
    all(isfinite, anchors) && all((0.0 .< anchors) .& (anchors .<= 1.0)) || throw(ArgumentError("anchor_pext must be finite and in (0, 1]"))

    iucn_levels = rand(rng, 1:5, Int(n_tips))
    tip_pext = anchors[iucn_levels]
    perm = sortperm(tip_pext)
    sorted_pext = tip_pext[perm]
    ranks = Float64.(1:Int(n_tips))
    X = hcat(ones(Int(n_tips)), ranks, ranks .^ 2, ranks .^ 3, ranks .^ 4)
    fitted = clamp.(X * (X \ sorted_pext), 0.0001, 0.9999)

    categories = String[]
    values = Float64[]
    for (name, level) in zip(("LC", "NT", "VU", "EN", "CR"), 1:5)
        for value in fitted[iucn_levels .== level]
            push!(categories, name)
            push!(values, value)
        end
    end
    return PextDistribution(categories, values)
end

function _edge2_column(df::AbstractDataFrame, wanted::AbstractString)
    matches = findall(name -> lowercase(String(name)) == lowercase(wanted), names(df))
    length(matches) == 1 || throw(ArgumentError("CSV must contain exactly one `$wanted` column"))
    return names(df)[only(matches)]
end

"""Load a pext distribution CSV with `RL.cat` and `pext` columns."""
function load_pext_distribution_csv(filepath::AbstractString)
    df = DataFrame(CSV.File(filepath))
    category_col = _edge2_column(df, "RL.cat")
    pext_col = _edge2_column(df, "pext")
    return PextDistribution(String.(df[!, category_col]), Float64.(df[!, pext_col]))
end

"""Sample extinction probabilities for IUCN-style category labels."""
function sample_pext(species_categories::AbstractVector{<:AbstractString}, distribution::PextDistribution; rng::AbstractRNG = Random.default_rng())
    pools = Dict(
        "LC" => _edge2_category_pool(distribution, ["LC"]),
        "NT" => _edge2_category_pool(distribution, ["NT"]),
        "VU" => _edge2_category_pool(distribution, ["VU"]),
        "EN" => _edge2_category_pool(distribution, ["EN"]),
        "CR" => _edge2_category_pool(distribution, ["CR"]),
        "thr" => _edge2_category_pool(distribution, ["VU", "EN", "CR"]),
        "not" => _edge2_category_pool(distribution, ["LC", "NT"]),
    )
    output = Vector{Float64}(undef, length(species_categories))
    @inbounds for (i, category0) in enumerate(species_categories)
        category = String(category0)
        if category == "EX"
            output[i] = 1.0
        elseif category in ("CR", "EW")
            output[i] = rand(rng, pools["CR"])
        elseif haskey(pools, category)
            output[i] = rand(rng, pools[category])
        elseif category == "thr" || category == "not"
            output[i] = rand(rng, pools[category])
        else
            output[i] = rand(rng, distribution.pext_values)
        end
    end
    return output
end

"""Reorder CSV species probabilities to `CompactTree.tip_labels` order."""
function reorder_pext_to_tree(tree::CompactTree, csv_species::AbstractVector{<:AbstractString}, csv_pext::AbstractVector{<:Real})
    length(csv_species) == length(csv_pext) || throw(ArgumentError("csv_species and csv_pext must have the same length"))
    lookup = Dict{String, Float64}()
    for (species, value) in zip(csv_species, csv_pext)
        haskey(lookup, String(species)) && throw(ArgumentError("duplicate species in pext data: $species"))
        lookup[String(species)] = Float64(value)
    end
    missing_species = setdiff(tree.tip_labels, collect(keys(lookup)))
    isempty(missing_species) || throw(ArgumentError("pext data is missing tree tips: $(join(missing_species, ", "))"))
    return [lookup[species] for species in tree.tip_labels]
end

"""Shift pext values by a calibrated category amount."""
function score_adjust(initial_scores::AbstractVector{<:Real}, distribution::PextDistribution, category_increase::Real)
    scores = Float64.(initial_scores)
    all(isfinite, scores) && all((0.0 .< scores) .& (scores .<= 1.0)) || throw(ArgumentError("initial_scores must be finite and in (0, 1]"))
    category_number = Dict("LC" => 1.0, "NT" => 2.0, "VU" => 3.0, "EN" => 4.0, "CR" => 5.0)
    calibrated_pext = sort(unique(distribution.pext_values))
    calibrated_category = Vector{Float64}(undef, length(calibrated_pext))
    for (i, value) in enumerate(calibrated_pext)
        indices = findall(x -> isapprox(x, value; atol = 1e-12), distribution.pext_values)
        counts = Dict{String, Int}()
        for category in distribution.categories[indices]
            counts[category] = get(counts, category, 0) + 1
        end
        best_category = isempty(counts) ? "VU" : first(sort!(collect(keys(counts)); by = category -> (-counts[category], category)))
        calibrated_category[i] = get(category_number, best_category, 3.0)
    end

    order = sortperm(scores)
    sorted_scores = scores[order]
    categories = similar(sorted_scores)
    for i in eachindex(sorted_scores)
        if sorted_scores[i] <= calibrated_pext[1]
            categories[i] = calibrated_category[1]
        elseif sorted_scores[i] >= calibrated_pext[end]
            categories[i] = calibrated_category[end]
        else
            index = searchsortedlast(calibrated_pext, sorted_scores[i])
            fraction = (sorted_scores[i] - calibrated_pext[index]) / (calibrated_pext[index + 1] - calibrated_pext[index])
            categories[i] = calibrated_category[index] + fraction * (calibrated_category[index + 1] - calibrated_category[index])
        end
    end
    shifted = clamp.(categories .+ Float64(category_increase), 0.0, 5.0)
    category_order = sortperm(calibrated_category)
    sorted_categories = calibrated_category[category_order]
    sorted_pext = calibrated_pext[category_order]
    adjusted = similar(shifted)
    for i in eachindex(shifted)
        if shifted[i] <= sorted_categories[1]
            adjusted[i] = sorted_pext[1]
        elseif shifted[i] >= sorted_categories[end]
            adjusted[i] = sorted_pext[end]
        else
            index = searchsortedlast(sorted_categories, shifted[i])
            fraction = (shifted[i] - sorted_categories[index]) / (sorted_categories[index + 1] - sorted_categories[index])
            adjusted[i] = sorted_pext[index] + fraction * (sorted_pext[index + 1] - sorted_pext[index])
        end
    end
    adjusted = clamp.(adjusted, 0.0001, 0.9999)
    @inbounds for i in eachindex(scores)
        scores[i] == 1.0 && (adjusted[findfirst(==(i), order)] = 1.0)
    end
    result = similar(adjusted)
    result[order] = adjusted
    return result
end
