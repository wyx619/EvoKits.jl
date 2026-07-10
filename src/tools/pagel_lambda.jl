Base.@kwdef struct PagelLambdaResult
    success::Bool = false
    lambda::Float64 = NaN
    loglik::Float64 = -Inf
    aic::Float64 = Inf
    sigma2::Float64 = NaN
    n::Int = 0
    converged::Bool = false
    iterations::Int = 0
    f_calls::Int = 0
    ci_lower::Float64 = NaN
    ci_upper::Float64 = NaN
    lr_stat::Float64 = NaN
    p_value::Float64 = NaN
    fit::Any = nothing
end

Base.@kwdef struct _PagelGLSProfile
    success::Bool = false
    loglik::Float64 = -Inf
    sigma2hat::Float64 = NaN
    message::String = ""
end

function _pagel_tip_crossproducts(tree::CompactTree, z::AbstractMatrix{<:Real}, edge_v::AbstractVector{<:Real})
    size(z, 1) == tree.ntips || throw(ArgumentError("z must have one row per tip"))
    length(edge_v) == tree.nedges || throw(ArgumentError("edge_v must have $(tree.nedges) entries"))
    _validate_binary_tree(tree)
    k = size(z, 2)
    any(x -> !isfinite(x), z) && return (success = false, cross = zeros(k, k), logdet = Inf, message = "non-finite data")

    precision = zeros(Float64, tree.nnodes)
    linear = zeros(Float64, tree.nnodes, k)
    quad = zeros(Float64, tree.nnodes, k, k)
    work_quad = zeros(Float64, k, k)
    logdet = 0.0

    tip_index = zeros(Int, tree.nnodes)
    for (i, node) in enumerate(tree.tip_ids)
        tip_index[node] = i
    end

    for node in tree.postorder_internal
        precision[node] = 0.0
        @views linear[node, :] .= 0.0
        @views quad[node, :, :] .= 0.0

        for edge in tree.first_child_edge[node]:tree.last_child_edge[node]
            child = Int(tree.child_of_edge[edge])
            v = Float64(edge_v[edge])
            (isfinite(v) && v >= 0.0) || return (success = false, cross = zeros(k, k), logdet = Inf, message = "invalid edge variance")

            child_p = 0.0
            child_s = zeros(Float64, k)
            child_q = zeros(Float64, k, k)
            child_logdet = 0.0

            if tree.is_tip[child]
                vv = max(v, 1e-12)
                invv = 1.0 / vv
                @views child_s .= invv .* z[tip_index[child], :]
                @views child_q .= invv .* (z[tip_index[child], :] * z[tip_index[child], :]')
                child_p = invv
                child_logdet = log(vv)
            else
                p = precision[child]
                vv = max(v, 1e-12)
                denom = 1.0 + vv * p
                scale = 1.0 / denom
                child_p = p * scale
                @views child_s .= scale .* linear[child, :]
                @views child_q .= quad[child, :, :] .- (vv * scale) .* (linear[child, :] * linear[child, :]')
                child_logdet = log(denom)
            end

            precision[node] += child_p
            @views linear[node, :] .+= child_s
            @views quad[node, :, :] .+= child_q
            logdet += child_logdet
        end
    end

    root = Int(tree.root)
    p = precision[root]
    p > 0.0 || return (success = false, cross = zeros(k, k), logdet = Inf, message = "root precision is zero")
    scale = 1.0 / p
    @views work_quad .= quad[root, :, :] .- scale .* (linear[root, :] * linear[root, :]')
    @views quad[root, :, :] .= work_quad

    return (success = isfinite(logdet), cross = Symmetric(Matrix(@view quad[root, :, :])), logdet = logdet, message = "")
end

function _pagel_profile_gls_from_cross(cross::AbstractMatrix{<:Real}, logdet::Real, nobs::Integer)
    size(cross, 1) == 2 || throw(ArgumentError("Pagel lambda cross matrix must contain y and intercept"))
    yy = Float64(cross[1, 1])
    xy = Float64(cross[2, 1])
    xx = Float64(cross[2, 2])
    xx > 0.0 || return _PagelGLSProfile(success = false, message = "non-positive intercept precision")
    beta0 = xy / xx
    rss = yy - xy * beta0
    (isfinite(rss) && rss > 0.0 && isfinite(logdet)) || return _PagelGLSProfile(success = false, message = "non-positive profiled residual sum of squares")
    sigma2hat = rss / Int(nobs)
    loglik = -0.5 * (Int(nobs) * (log(2pi * sigma2hat) + 1.0) + Float64(logdet))
    return _PagelGLSProfile(success = isfinite(loglik), loglik = loglik, sigma2hat = sigma2hat)
end

function _pagel_lambda_edge_variances(tree::CompactTree, lambda::Real)
    lam = Float64(lambda)
    0.0 <= lam <= 1.0 || throw(ArgumentError("lambda must be within [0, 1]"))
    edge_v = lam .* Float64.(tree.edge_length)
    @inbounds for edge in 1:tree.nedges
        child = Int(tree.child_of_edge[edge])
        if tree.is_tip[child]
            edge_v[edge] += (1.0 - lam) * tree.dist_from_root[child]
        end
    end
    return edge_v
end

function _pagel_lambda_fit(tree::CompactTree, y::AbstractVector{<:Real}; max_iterations::Integer = 200, rel_tol::Float64 = 1e-8)
    yy = _validate_univariate_trait(tree, y)
    _validate_binary_tree(tree)
    _validate_ultrametric_tree(tree)
    z = hcat(Float64.(yy), ones(Float64, tree.ntips))
    eps_lambda = 1e-8

    objective = function (lambda)
        lam = clamp(Float64(lambda), eps_lambda, 1.0 - eps_lambda)
        edge_v = _pagel_lambda_edge_variances(tree, lam)
        cross = _pagel_tip_crossproducts(tree, z, edge_v)
        cross.success || return Inf
        prof = _pagel_profile_gls_from_cross(cross.cross, cross.logdet, tree.ntips)
        return prof.success ? -prof.loglik : Inf
    end

    result = Optim.optimize(objective, eps_lambda, 1.0 - eps_lambda, Optim.Brent(); iterations = max_iterations, rel_tol = rel_tol)
    lamhat = Optim.minimizer(result)
    edge_v = _pagel_lambda_edge_variances(tree, lamhat)
    cross = _pagel_tip_crossproducts(tree, z, edge_v)
    cross.success || return PagelLambdaResult(success = false, n = tree.ntips, converged = Optim.converged(result), iterations = Optim.iterations(result), f_calls = Optim.f_calls(result))
    prof = _pagel_profile_gls_from_cross(cross.cross, cross.logdet, tree.ntips)
    prof.success || return PagelLambdaResult(success = false, n = tree.ntips, converged = Optim.converged(result), iterations = Optim.iterations(result), f_calls = Optim.f_calls(result))

    return PagelLambdaResult(
        success = true,
        lambda = lamhat,
        loglik = prof.loglik,
        aic = compute_aic(prof.loglik, 2),
        sigma2 = prof.sigma2hat,
        n = tree.ntips,
        converged = Optim.converged(result),
        iterations = Optim.iterations(result),
        f_calls = Optim.f_calls(result),
        fit = (model = :lambda, lambda = lamhat, sigma2 = prof.sigma2hat, loglik = prof.loglik, aic = compute_aic(prof.loglik, 2)),
    )
end

@inline function _pagel_lambda_lr_pvalue(lr_stat::Real)
    x = Float64(lr_stat)
    x <= 0.0 && return 1.0
    return ccdf(Chisq(1), x)
end

function _pagel_lambda_profile(tree::CompactTree, y::AbstractVector{<:Real}; npoints::Integer = 101)
    npoints >= 3 || throw(ArgumentError("npoints must be at least 3"))
    yy = _validate_univariate_trait(tree, y)
    z = hcat(Float64.(yy), ones(Float64, tree.ntips))
    eps_lambda = 1e-8
    lambdas = collect(range(eps_lambda, 1.0 - eps_lambda, length = Int(npoints)))
    logliks = fill(-Inf, length(lambdas))

    for i in eachindex(lambdas)
        edge_v = _pagel_lambda_edge_variances(tree, lambdas[i])
        cross = _pagel_tip_crossproducts(tree, z, edge_v)
        cross.success || continue
        prof = _pagel_profile_gls_from_cross(cross.cross, cross.logdet, tree.ntips)
        prof.success || continue
        logliks[i] = prof.loglik
    end
    return lambdas, logliks
end

function _pagel_lambda_ci_from_profile(lambdas::AbstractVector{<:Real}, logliks::AbstractVector{<:Real}, mle::Real)
    llmax = maximum(logliks)
    isfinite(llmax) || return (NaN, NaN)
    threshold = llmax - 1.92
    mle_idx = findmin(abs.(Float64.(lambdas) .- Float64(mle)))[2]

    lower = 0.0
    upper = 1.0
    for i in 1:(mle_idx - 1)
        if logliks[i] < threshold && logliks[i + 1] >= threshold
            x1, x2 = Float64(lambdas[i]), Float64(lambdas[i + 1])
            y1, y2 = Float64(logliks[i]), Float64(logliks[i + 1])
            lower = x1 + (x2 - x1) * (threshold - y1) / (y2 - y1)
            break
        end
    end
    for i in mle_idx:(length(lambdas) - 1)
        if logliks[i] >= threshold && logliks[i + 1] < threshold
            x1, x2 = Float64(lambdas[i]), Float64(lambdas[i + 1])
            y1, y2 = Float64(logliks[i]), Float64(logliks[i + 1])
            upper = x1 + (x2 - x1) * (threshold - y1) / (y2 - y1)
            break
        end
    end
    return lower, upper
end

function pagel_lambda(tree::CompactTree, y::AbstractVector{<:Real}; root_model::Symbol = :random, max_iterations::Integer = 200, rel_tol::Float64 = 1e-8, profile_points::Integer = 101)
    root_model in (:random, :fixed) || throw(ArgumentError("root_model must be :random or :fixed"))
    fit = _pagel_lambda_fit(tree, y; max_iterations = max_iterations, rel_tol = rel_tol)
    !fit.success || isfinite(fit.lambda) || return fit

    lambdas, logliks = _pagel_lambda_profile(tree, y; npoints = profile_points)
    ci_lower, ci_upper = _pagel_lambda_ci_from_profile(lambdas, logliks, fit.lambda)

    edge_v0 = _pagel_lambda_edge_variances(tree, 1e-8)
    z = hcat(Float64.(_validate_univariate_trait(tree, y)), ones(Float64, tree.ntips))
    cross0 = _pagel_tip_crossproducts(tree, z, edge_v0)
    ll0 = if cross0.success
        prof0 = _pagel_profile_gls_from_cross(cross0.cross, cross0.logdet, tree.ntips)
        prof0.success ? prof0.loglik : -Inf
    else
        -Inf
    end
    lr_stat = isfinite(ll0) ? max(0.0, 2.0 * (fit.loglik - ll0)) : NaN
    p_value = isfinite(lr_stat) ? _pagel_lambda_lr_pvalue(lr_stat) : NaN

    return PagelLambdaResult(
        success = fit.success,
        lambda = fit.lambda,
        loglik = fit.loglik,
        aic = fit.aic,
        sigma2 = fit.sigma2,
        n = fit.n,
        converged = fit.converged,
        iterations = fit.iterations,
        f_calls = fit.f_calls,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        lr_stat = lr_stat,
        p_value = p_value,
        fit = fit.fit,
    )
end

function pagel_lambda(tree::CompactTree, data::AbstractDataFrame, trait_col; taxon_col::Union{Nothing, Symbol, AbstractString} = nothing, kwargs...)
    df = DataFrames.DataFrame(data)
    trait_sym = Symbol(trait_col)
    hasproperty(df, trait_sym) || throw(ArgumentError("data does not contain trait column `$trait_sym`"))
    if taxon_col === nothing
        DataFrames.nrow(df) == tree.ntips || throw(ArgumentError("data must have one row per tree tip when taxon_col is omitted"))
        y = Float64.(df[!, trait_sym])
        return pagel_lambda(tree, y; kwargs...)
    end
    aligned = align_traits_to_tree(tree, df; taxon_col = taxon_col, trait_cols = trait_sym)
    return pagel_lambda(tree, aligned; kwargs...)
end

function summary_pagel_lambda(res::PagelLambdaResult)
    return (
        success = res.success,
        lambda = res.lambda,
        logLik = res.loglik,
        AIC = res.aic,
        sigma2 = res.sigma2,
        n = res.n,
        converged = res.converged,
        iterations = res.iterations,
        f_calls = res.f_calls,
        ci_lower = res.ci_lower,
        ci_upper = res.ci_upper,
        lr_stat = res.lr_stat,
        p_value = res.p_value,
    )
end
