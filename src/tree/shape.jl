Base.@kwdef struct LineageThroughTimeResult
    event_times::Vector{Float64} = Float64[]
    lineages_after_event::Vector{Int} = Int[]
    event_nodes::Vector{Int32} = Int32[]
    query_times::Vector{Float64} = Float64[]
    lineages_at_query::Vector{Int} = Int[]
end

function tree_imbalance(tree::CompactTree; metric::Symbol = :sackin, normalized::Bool = false)
    counts = _tools_descendant_tip_counts(tree)
    if metric === :sackin
        depths = distances_from_root(tree; weighted = false)
        return Float64(sum(depths[Int(tip)] for tip in tree.tip_ids))
    elseif metric === :blum
        return sum(log(counts[Int(node)] - 1) for node in tree.postorder_internal if counts[Int(node)] > 1)
    elseif metric === :colless || metric === :colless_normalized
        value = 0.0
        for node32 in tree.postorder_internal
            child_counts = [counts[Int(child)] for child in tree.children[Int(node32)]]
            for i in 1:length(child_counts), j in (i + 1):length(child_counts)
                value += abs(child_counts[i] - child_counts[j])
            end
        end
        if normalized || metric === :colless_normalized
            denominator = 0.5 * (tree.ntips - 1) * (tree.ntips - 2)
            value = denominator > 0.0 ? value / denominator : 0.0
        end
        return value
    else
        throw(ArgumentError("metric must be :sackin, :colless, :colless_normalized, or :blum"))
    end
end

function node_depths(tree::CompactTree; weighted::Bool = true)
    return distances_from_root(tree; weighted = weighted)
end

function gamma_statistic(tree::CompactTree; atol::Float64 = 1e-8)
    tip_depths = tree.dist_from_root[tree.tip_ids]
    maximum(abs.(tip_depths .- first(tip_depths))) <= atol ||
        throw(ArgumentError("gamma_statistic requires an ultrametric tree"))
    nnodes = length(tree.postorder_internal)
    nnodes >= 2 || return NaN
    crown_height = maximum(tip_depths)
    ages = crown_height .- tree.dist_from_root[tree.postorder_internal]
    order = sortperm(ages; rev = true)
    event_nodes = tree.postorder_internal[order]
    event_ages = ages[order]
    lineages = Vector{Int}(undef, nnodes)
    for i in 1:nnodes
        lineages[i] = (i == 1 ? 1 : lineages[i - 1]) + length(tree.children[Int(event_nodes[i])]) - 1
    end
    event_times = zeros(Float64, nnodes)
    for i in 1:nnodes
        gap = event_ages[i] - (i == nnodes ? 0.0 : event_ages[i + 1])
        event_times[i] = (i == 1 ? 0.0 : event_times[i - 1]) + lineages[i] * gap
    end
    total = event_times[end]
    total > 0.0 || return NaN
    return ((sum(@view(event_times[1:(end - 1)])) / (nnodes - 1)) - 0.5 * total) /
        (total / sqrt(12.0 * (nnodes - 1)))
end

function lineages_through_time(tree::CompactTree; times = nothing, atol::Float64 = 1e-8)
    tip_depths = tree.dist_from_root[tree.tip_ids]
    maximum(abs.(tip_depths .- first(tip_depths))) <= atol ||
        throw(ArgumentError("lineages_through_time requires an ultrametric tree"))
    height = maximum(tip_depths)
    nodes = copy(tree.postorder_internal)
    node_times = [tree.dist_from_root[Int(node)] for node in nodes]
    order = sortperm(node_times)
    ordered_nodes = nodes[order]
    ordered_times = node_times[order]
    lineages = Vector{Int}(undef, length(ordered_nodes))
    current = 1
    for i in eachindex(ordered_nodes)
        current += length(tree.children[Int(ordered_nodes[i])]) - 1
        lineages[i] = current
    end
    query_times = times === nothing ? Float64[] : Float64.(collect(times))
    query_values = Int[]
    for time in query_times
        0.0 <= time <= height || throw(ArgumentError("query times must lie within 0 and tree height"))
        value = 1 + sum(length(tree.children[Int(node)]) - 1 for (node, node_time) in zip(nodes, node_times) if node_time <= time)
        push!(query_values, value)
    end
    return LineageThroughTimeResult(
        event_times = ordered_times,
        lineages_after_event = lineages,
        event_nodes = ordered_nodes,
        query_times = query_times,
        lineages_at_query = query_values,
    )
end
