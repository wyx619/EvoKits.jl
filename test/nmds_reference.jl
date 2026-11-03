# Independent weighted min-max characterization of isotonic regression.
# This deliberately uses BigFloat interval means, not PAVA or production code.
function nmds_reference_fit(dissimilarities, distances, ties)
    order = sortperm(eachindex(distances); by = i -> ties === :weak ?
        (dissimilarities[i], distances[i], i) : (dissimilarities[i], 0.0, i))
    groups = Vector{Int}[]
    for i in order
        if ties === :strong && !isempty(groups) && dissimilarities[first(last(groups))] == dissimilarities[i]
            push!(last(groups), i)
        else
            push!(groups, [i])
        end
    end
    sums = [sum(BigFloat(distances[i]) for i in group) for group in groups]
    weights = length.(groups)
    fitted = similar(distances)
    for i in eachindex(groups)
        value = maximum(minimum(sum(sums[a:b]) / sum(weights[a:b])
            for b in i:length(groups)) for a in 1:i)
        fitted[groups[i]] .= Float64(value)
    end
    return fitted
end

function nmds_reference_distances(points)
    n = size(points, 1)
    return [Float64(sqrt(sum(abs2, BigFloat.(points[i, :]) .- BigFloat.(points[j, :]))))
        for j in 2:n for i in 1:(j - 1)]
end

function nmds_reference_stress(observed, points, ties)
    distances = nmds_reference_distances(points)
    fitted = nmds_reference_fit(observed, distances, ties)
    numerator = sum(abs2, BigFloat.(distances) .- BigFloat.(fitted))
    denominator = sum(abs2, BigFloat.(distances))
    return (; distances, fitted, stress = Float64(sqrt(numerator / denominator)))
end

function nmds_reference_step(observed, points, ties)
    n = size(points, 1)
    fit = nmds_reference_stress(observed, points, ties)
    b = zeros(BigFloat, n, n)
    index = 0
    for j in 2:n, i in 1:(j - 1)
        index += 1
        ratio = fit.distances[index] == 0 ? BigFloat(0) :
            BigFloat(fit.fitted[index]) / BigFloat(fit.distances[index])
        b[i, j] = b[j, i] = -ratio
        b[i, i] += ratio
        b[j, j] += ratio
    end
    updated = b * BigFloat.(points) / n
    for axis in axes(updated, 2)
        updated[:, axis] .-= sum(updated[:, axis]) / n
    end
    updated ./= sqrt(sum(abs2, updated) / n)
    return Float64.(updated)
end

function nmds_test_matrix(packed, n)
    matrix = zeros(n, n)
    index = 0
    for j in 2:n, i in 1:(j - 1)
        index += 1
        matrix[i, j] = matrix[j, i] = packed[index]
    end
    return matrix
end
