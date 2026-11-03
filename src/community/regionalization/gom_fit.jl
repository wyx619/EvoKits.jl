# Scalar-shape NEF objective and EM updates used by maptpx. Responsibilities
# are accumulated one observed entry at a time, without an nnz × k array.
function _gom_counts(community::AbstractMatrix{<:Real})
    n, p = size(community)
    n > 0 && p > 0 || throw(ArgumentError("fitgom requires sites and species"))
    _validate_community_values(community)
    rows, columns, values = Int[], Int[], Float64[]
    row_totals = zeros(n)
    total = 0.0
    _community_each_positive(community) do i, j, value
        x = Float64(value)
        isfinite(x) && x > 0 ||
            throw(ArgumentError("positive counts must be representable as finite positive Float64 values"))
        total += x
        row_totals[i] += x
        isfinite(total) || throw(ArgumentError("count totals exceed Float64 range; rescale counts"))
        push!(rows, i)
        push!(columns, j)
        push!(values, x)
    end
    all(>(0), row_totals) ||
        throw(ArgumentError("fitgom does not assign empty sites; select nonempty sites and their labels explicitly"))
    return (; rows, columns, values, row_totals, total, n, p)
end

# Return log(q_ij), writing the k responsibilities into a reusable vector.
function _gom_edge!(work, theta, omega, i, j)
    q = 0.0
    @inbounds for r in eachindex(work)
        work[r] = theta[j, r] * omega[i, r]
        q += work[r]
    end
    if q >= floatmin(Float64)
        @inbounds for r in eachindex(work)
            work[r] /= q
        end
        return log(q)
    end
    # Strictly positive factors can still have products below Float64 range.
    largest = -Inf
    @inbounds for r in eachindex(work)
        work[r] = log(theta[j, r]) + log(omega[i, r])
        largest = max(largest, work[r])
    end
    normalizer = 0.0
    @inbounds for r in eachindex(work)
        work[r] = exp(work[r] - largest)
        normalizer += work[r]
    end
    @inbounds for r in eachindex(work)
        work[r] /= normalizer
    end
    return largest + log(normalizer)
end

function _gom_score(counts, theta, omega, shape, work)
    likelihood = 0.0
    @inbounds for e in eachindex(counts.values)
        likelihood += counts.values[e] *
            _gom_edge!(work, theta, omega, counts.rows[e], counts.columns[e])
    end
    objective = likelihood + shape * sum(log, theta) + sum(log, omega) / size(omega, 2)
    return likelihood, objective
end

function _gom_em!(next_theta, next_omega, counts, theta, omega, shape, work)
    k = size(theta, 2)
    fill!(next_theta, shape)
    fill!(next_omega, 1 / k)
    @inbounds for e in eachindex(counts.values)
        i, j, x = counts.rows[e], counts.columns[e], counts.values[e]
        _gom_edge!(work, theta, omega, i, j)
        for r in 1:k
            z = x * work[r]
            next_theta[j, r] += z
            next_omega[i, r] += z
        end
    end
    for r in 1:k
        column = @view next_theta[:, r]
        column ./= sum(column)
    end
    for i in 1:counts.n
        row = @view next_omega[i, :]
        row ./= sum(row)
    end
    return nothing
end

_gom_difference(a, b) = maximum(abs(a[i] - b[i]) for i in eachindex(a, b))

function _gom_initial(input, n, k, byrow, name)
    input isa AbstractMatrix{<:Real} && size(input) == (n, k) ||
        throw(ArgumentError("$name must be a $n × $k matrix of positive real values"))
    result = Matrix{Float64}(undef, n, k)
    for i in eachindex(result, input)
        x = input[i]
        value = Float64(x)
        isfinite(x) && x > 0 && isfinite(value) && value > 0 ||
            throw(ArgumentError("$name must contain finite positive Float64 values"))
        result[i] = value
    end
    groups = byrow ? eachrow(result) : eachcol(result)
    for group in groups
        group ./= maximum(group)
        group ./= sum(group)
    end
    all(>(0), result) || throw(ArgumentError("normalizing $name underflows; use less extreme initial values"))
    return result
end

function _gom_random_start(counts, k, rng)
    theta = 1 .+ rand(rng, counts.p, k)
    for column in eachcol(theta)
        column .*= 0.1 / sum(column)
    end
    anchors = zeros(Int, counts.n)
    for (r, i) in enumerate(randperm(rng, counts.n)[1:k])
        anchors[i] = r
    end
    @inbounds for e in eachindex(counts.values)
        i, j = counts.rows[e], counts.columns[e]
        r = anchors[i]
        r == 0 && continue
        theta[j, r] += 0.9 * (counts.values[e] / counts.row_totals[i])
    end
    for column in eachcol(theta)
        column ./= sum(column)
    end
    omega = 1 .+ rand(rng, counts.n, k)
    for row in eachrow(omega)
        row ./= sum(row)
    end
    return theta, omega
end

# Work in log ratios with a reference category, so extrapolation preserves the
# simplexes after softmax. The secant proposal is accepted only if it improves
# on both ordinary EM steps. This omits maptpx's conditional SQP omega solver.
function _gom_nef!(yt, yw, theta, omega)
    for r in axes(theta, 2), j in axes(theta, 1)
        yt[j, r] = log(theta[j, r]) - log(theta[1, r])
    end
    for r in axes(omega, 2), i in axes(omega, 1)
        yw[i, r] = log(omega[i, r]) - log(omega[i, 1])
    end
end

function _gom_simplex!(matrix, byrow)
    groups = byrow ? eachrow(matrix) : eachcol(matrix)
    for group in groups
        largest = maximum(group)
        isfinite(largest) || return false
        group .= exp.(group .- largest)
        group ./= sum(group)
    end
    return all(>(0), matrix)
end

function _gom_fit(counts, theta, omega, shape, maxiter, tol, accelerate)
    first_theta, next_theta, candidate_theta = similar(theta), similar(theta), similar(theta)
    first_omega, next_omega, candidate_omega = similar(omega), similar(omega), similar(omega)
    yt = [similar(theta) for _ in 1:3]
    yw = [similar(omega) for _ in 1:3]
    work = zeros(size(theta, 2))
    likelihood, objective = _gom_score(counts, theta, omega, shape, work)
    isfinite(objective) || throw(ArgumentError("initial objective exceeds Float64 range; rescale counts or priors"))
    history = [objective]
    reason, iterations = :maxiter, 0
    for iteration in 1:maxiter
        _gom_em!(first_theta, first_omega, counts, theta, omega, shape, work)
        _gom_em!(next_theta, next_omega, counts, first_theta, first_omega, shape, work)
        next_likelihood, next_objective = _gom_score(counts, next_theta, next_omega, shape, work)
        residual = max(_gom_difference(first_theta, next_theta), _gom_difference(first_omega, next_omega))
        if !isfinite(next_objective) || next_objective < objective
            # Preserve monotonic history even at the floating-point floor.
            reason = :stalled
            break
        end
        done = residual <= tol && next_objective - objective <= tol * (1 + abs(objective))
        if accelerate && !done
            _gom_nef!(yt[1], yw[1], theta, omega)
            _gom_nef!(yt[2], yw[2], first_theta, first_omega)
            _gom_nef!(yt[3], yw[3], next_theta, next_omega)
            rr, rv = 0.0, 0.0
            for states in (yt, yw), i in eachindex(states[1])
                r = states[2][i] - states[1][i]
                v = states[3][i] - states[2][i]
                rr += r * r
                rv += r * v
            end
            factor = rv / (rr - rv)
            if isfinite(factor) && 0 < factor <= 1000
                @. candidate_theta = yt[3] + factor * (yt[3] - yt[2])
                @. candidate_omega = yw[3] + factor * (yw[3] - yw[2])
                if _gom_simplex!(candidate_theta, false) && _gom_simplex!(candidate_omega, true)
                    candidate_likelihood, candidate_objective =
                        _gom_score(counts, candidate_theta, candidate_omega, shape, work)
                    if isfinite(candidate_objective) && candidate_objective > next_objective
                        next_theta .= candidate_theta
                        next_omega .= candidate_omega
                        next_likelihood, next_objective = candidate_likelihood, candidate_objective
                    end
                end
            end
        end
        theta, next_theta = next_theta, theta
        omega, next_omega = next_omega, omega
        likelihood, objective = next_likelihood, next_objective
        iterations = iteration
        push!(history, objective)
        if done
            # EM changes can grow near a symmetric stationary point. Check
            # the accepted factors before stopping on the preliminary change.
            _gom_em!(first_theta, first_omega, counts, theta, omega, shape, work)
            if max(_gom_difference(theta, first_theta), _gom_difference(omega, first_omega)) <= tol
                reason = :tolerance
                break
            end
        end
    end
    _gom_em!(first_theta, first_omega, counts, theta, omega, shape, work)
    residual = max(_gom_difference(theta, first_theta), _gom_difference(omega, first_omega))
    converged = residual <= tol
    converged && (reason = :tolerance)
    return (; theta, omega, likelihood, objective, iterations, converged, reason, residual, history)
end
