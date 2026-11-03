# Reuse the independent dense sampler used to validate the existing sparse
# null engine; it does not call any production randomizer or SES accumulator.
if !isdefined(@__MODULE__, :ses_dense_draw)
    include("ses_reference.jl")
end

function beta_ses_reference_values(presence, tree, family, component)
    # Compute descendant occupancy, independently of the production ancestor
    # walks. BigFloat sums separate common and unique branches directly.
    descendants = [Int[] for _ in 1:tree.nnodes]
    tip_positions = Dict(node => i for (i, node) in enumerate(tree.tip_ids))
    for node in tree.postorder
        if tree.is_tip[node]
            push!(descendants[node], tip_positions[node])
        else
            for child in tree.children[node]
                append!(descendants[node], descendants[child])
            end
        end
    end
    n = size(presence, 1)
    active = [any(presence[i, tip] for tip in descendants[tree.child_of_edge[e]])
        for i in 1:n, e in 1:tree.nedges]
    result = fill(NaN, n, n)
    for j in 1:n, i in j:n
        a, b, c = BigFloat(0), BigFloat(0), BigFloat(0)
        for e in 1:tree.nedges
            left, right = active[i, e], active[j, e]
            if left && right
                a += BigFloat(tree.edge_length[e])
            elseif left
                b += BigFloat(tree.edge_length[e])
            elseif right
                c += BigFloat(tree.edge_length[e])
            end
        end
        total = family === :sorensen ? (b + c) / (2a + b + c) : (b + c) / (a + b + c)
        turnover = family === :sorensen ? min(b, c) / (a + min(b, c)) :
            2min(b, c) / (a + 2min(b, c))
        value = component === :total ? total : component === :turnover ? turnover : total - turnover
        result[i, j] = result[j, i] = Float64(value)
    end
    return result
end

function beta_ses_reference(community, tree, model, family, component;
    runs = 13, iterations = 15, attempts = 35, seed = 42, pairs = nothing,
)
    n = size(community, 1)
    selected = pairs === nothing ? [(i, j) for j in 1:n for i in 1:n] : pairs
    presence = community .> 0
    rng = StableRNG(seed)
    exact_obs = beta_ses_reference_values(presence, tree, family, component)
    public_obs = getproperty(phylobeta(presence, tree; index_family = family), component)
    draws = Matrix{Float64}(undef, length(selected), runs)
    rank_draws = similar(draws)
    for draw in 1:runs
        random = ses_dense_draw(presence, model, iterations, attempts, rng)
        independent = beta_ses_reference_values(random, tree, family, component)
        public = getproperty(phylobeta(random, tree; index_family = family), component)
        for (k, (i, j)) in enumerate(selected)
            draws[k, draw] = independent[i, j]
            rank_draws[k, draw] = public[i, j]
        end
    end
    means, sds, z, p = (fill(NaN, length(selected)) for _ in 1:4)
    counts = zeros(Int, length(selected))
    for (k, (i, j)) in enumerate(selected)
        finite = filter(isfinite, draws[k, :])
        counts[k] = length(finite)
        if !isempty(finite)
            bigdraws = BigFloat.(finite)
            mu = sum(bigdraws) / length(bigdraws)
            means[k] = Float64(mu)
            if length(finite) > 1
                sds[k] = Float64(sqrt(sum(abs2, bigdraws .- mu) / (length(bigdraws) - 1)))
                sds[k] > 0 && (z[k] = (exact_obs[i, j] - means[k]) / sds[k])
            end
            if isfinite(public_obs[i, j])
                rank_values = filter(isfinite, rank_draws[k, :])
                p[k] = (1 + count(<(public_obs[i, j]), rank_values) +
                    0.5count(==(public_obs[i, j]), rank_values)) / (length(rank_values) + 1)
            end
        end
    end
    return (; observed = [exact_obs[i, j] for (i, j) in selected], null_mean = means,
        null_sd = sds, z_score = z, p_value = p, valid_replicates = counts,
        next_random = rand(rng, UInt64))
end
