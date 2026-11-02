using Statistics

# Independent dense null-model oracle. Full permutation calls and checkerboard
# trials reproduce the pre-refactor draws, without using production samplers.
function ses_dense_draw(presence, model, iterations, attempts, rng)
    n, p = size(presence)
    result = fill(false, n, p)
    if model === :taxa_labels
        result .= presence[:, randperm(rng, p)]
    elseif model in (:richness, :sample_pool, :phylogeny_pool)
        for i in 1:n
            richness = count(presence[i, :])
            richness == 0 && continue
            result[i, randperm(rng, p)[1:richness]] .= true
        end
        model === :phylogeny_pool && (result = result[:, randperm(rng, p)])
    elseif model === :frequency
        for j in 1:p
            frequency = count(presence[:, j])
            frequency == 0 && continue
            result[randperm(rng, n)[1:frequency], j] .= true
        end
    else
        result .= presence
        if n >= 2 && p >= 2
            function pick_pair(size)
                a, b = rand(rng, 1:size), rand(rng, 1:size-1)
                return a, b >= a ? b + 1 : b
            end
            for _ in 1:iterations
                swapped = false
                for _ in 1:(model === :trialswap ? 1 : attempts)
                    i, j = pick_pair(n)
                    k, l = pick_pair(p)
                    a, b, c, d = result[i, k], result[i, l], result[j, k], result[j, l]
                    if a == d && b == c && a != b
                        result[i, k], result[i, l] = b, a
                        result[j, k], result[j, l] = d, c
                        swapped = true
                        break
                    end
                end
                !swapped && model === :independentswap && break
            end
        end
    end
    return result
end

function ses_metric_reference(presence, abundance, tree, metric; weighted = false, include_root = true)
    if metric === :pd
        reference = community_reference(presence, tree)
        return include_root ? reference.pd : reference.pd_no_root
    end
    distances = pairwise_distance(tree)
    values = fill(NaN, size(presence, 1))
    for i in axes(presence, 1)
        present = findall(presence[i, :])
        isempty(present) && continue
        d = distances[present, present]
        weights = weighted ? abundance[i, present] : ones(length(present))
        total = sum(weights)
        total > 0 || continue
        if metric === :mpd
            if weighted
                values[i] = sum(d .* (weights * weights')) / total^2
            elseif length(present) > 1
                values[i] = sum(d) / (length(present) * (length(present) - 1))
            end
        elseif length(present) > 1
            nearest = [minimum(d[j, k] for k in axes(d, 2) if j != k) for j in axes(d, 1)]
            values[i] = sum(nearest .* weights) / total
        end
    end
    return values
end

function ses_summary_reference(community, tree, metric, model; runs = 11, iterations = 20,
    attempts = 40, weighted = false, include_root = true, seed = 123)
    rng = MersenneTwister(seed)
    presence = community .> 0
    obs = ses_metric_reference(presence, community, tree, metric; weighted, include_root)
    draws = Matrix{Float64}(undef, size(community, 1), runs)
    for run in 1:runs
        random = ses_dense_draw(presence, model, iterations, attempts, rng)
        draws[:, run] = ses_metric_reference(random, community, tree, metric; weighted, include_root)
    end
    means, sds, ranks, z = (fill(NaN, length(obs)) for _ in 1:4)
    for i in eachindex(obs)
        finite = filter(isfinite, draws[i, :])
        # Center the independent two-pass variance at an actual draw. Summing
        # identical decimals first can create a spurious nonzero SD through
        # rounding of their mean (e.g. eleven copies of 4.4).
        offsets = isempty(finite) ? Float64[] : finite .- first(finite)
        offset_mean = isempty(offsets) ? NaN : mean(offsets)
        means[i] = isempty(finite) ? NaN : first(finite) + offset_mean
        sds[i] = length(finite) > 1 ? sqrt(sum(abs2, offsets .- offset_mean) / (length(finite) - 1)) : NaN
        isfinite(sds[i]) && sds[i] > 0 && (z[i] = (obs[i] - means[i]) / sds[i])
        ranks[i] = 1 + count(x -> isfinite(obs[i]) && x < obs[i], finite) +
            0.5count(x -> isfinite(obs[i]) && x == obs[i], finite)
    end
    return (; obs, means, sds, ranks, z, p = ranks ./ (runs + 1), draws)
end
