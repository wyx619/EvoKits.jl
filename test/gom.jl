using Test, EvoKits, SparseArrays, Random, StableRNGs, LinearAlgebra, CSV, DataFrames

# Dense BigFloat oracle written from the model, including unobserved cells.
function gom_test_em(x, theta, omega, shape)
    x, theta, omega = BigFloat.(x), BigFloat.(theta), BigFloat.(omega)
    n, p = size(x)
    k = size(theta, 2)
    q = omega * transpose(theta)
    z = [x[i,j] * omega[i,r] * theta[j,r] / q[i,j] for i in 1:n, j in 1:p, r in 1:k]
    next_theta = [sum(z[:,j,r]) + shape for j in 1:p, r in 1:k]
    next_omega = [sum(z[i,:,r]) + inv(BigFloat(k)) for i in 1:n, r in 1:k]
    next_theta ./= sum(next_theta; dims=1)
    next_omega ./= sum(next_omega; dims=2)
    likelihood = sum(x .* log.(q))
    objective = likelihood + BigFloat(shape) * sum(log, theta) + sum(log, omega) / k
    return Float64.(next_theta), Float64.(next_omega), Float64(likelihood), Float64(objective)
end

@testset "GoM: independent objective and EM pseudocounts" begin
    rng = StableRNG(206)
    for k in 1:3, shape in (0.125, 1.7), sparse_input in (false, true)
        x = rand(rng, 0:9, 4, 5)
        x[:,1] .+= 1
        x[:,5] .= 0
        theta = rand(rng, 5, k) .+ 0.1
        theta ./= sum(theta; dims=1)
        omega = rand(rng, 4, k) .+ 0.1
        omega ./= sum(omega; dims=2)
        expected_t, expected_w, likelihood, objective = gom_test_em(x, theta, omega, shape)
        counts = EvoKits._gom_counts(sparse_input ? sparse(x) : x)
        actual_t, actual_w = similar(theta), similar(omega)
        EvoKits._gom_em!(actual_t, actual_w, counts, theta, omega, shape, zeros(k))
        @test actual_t ≈ expected_t atol=1e-14
        @test actual_w ≈ expected_w atol=1e-14
        actual_l, actual_objective = EvoKits._gom_score(counts, theta, omega, shape, zeros(k))
        @test actual_l ≈ likelihood atol=1e-12
        @test actual_objective ≈ objective atol=1e-12
        @test all(>(0), actual_t) && all(>(0), actual_w)
        @test sum(actual_t; dims=1) ≈ ones(1,k)
        @test sum(actual_w; dims=2) ≈ ones(4,1)
    end
    # EM log-sum-exp path: every ordinary product for this observation underflows.
    theta = [1e-200 1.0; 1.0 1e-200]
    omega = [1.0 1e-200; 1e-200 1.0]
    work = zeros(2)
    logq = EvoKits._gom_edge!(work, theta, omega, 1, 1)
    @test logq ≈ log(2.0) - 200log(10.0)
    @test work ≈ [0.5, 0.5]
    tiny_t = [1e-200 1e-200; 1.0 1.0]
    tiny_w = [1e-200 1e-200; 1.0 1.0]
    @test EvoKits._gom_edge!(work, tiny_t, tiny_w, 1, 1) ≈ log(2.0) - 400log(10.0)
    @test work ≈ [0.5, 0.5]
end

@testset "GoM: fit, sparse parity, diagnostics and reproducibility" begin
    x = [20 10 0 0; 15 25 0 0; 0 0 20 10; 0 0 15 25; 10 5 5 10; 5 10 10 5]
    samples = ["west_a", "west_b", "east_a", "east_b", "mixed_a", "mixed_b"]
    species = ["A", "B", "C", "D"]
    data = CommunityData(sparse(x); sample_labels=samples, species_labels=species)
    before = copy(data.matrix)
    for accelerate in (false, true)
        fit = fitgom(data; k=2, nstart=3, rng=StableRNG(42), accelerate)
        dense = fitgom(x; k=2, nstart=3, rng=StableRNG(42), accelerate,
            sample_labels=samples, species_labels=species)
        @test fit isa GoMResult
        @test fit.omega == dense.omega && fit.theta == dense.theta
        @test fit.history == dense.history && fit.start_objectives == dense.start_objectives
        @test fit.shape == 1/8 && fit.k == 2
        @test fit.maxiter == 1000 && fit.tol == 1e-8 && fit.accelerate == accelerate
        @test fit.sample_labels == samples && fit.species_labels == species
        @test fit.sample_labels !== data.sample_labels && fit.species_labels !== data.species_labels
        @test fit.region_labels == ["region_1", "region_2"]
        @test fit.converged && fit.reason === :tolerance && fit.residual <= 1e-8
        @test length(fit.history) == fit.iterations + 1
        @test all(>=(0), diff(fit.history)) && last(fit.history) == fit.objective
        @test fit.objective == maximum(fit.start_objectives) == fit.start_objectives[fit.best_start]
        @test fit.start_iterations[fit.best_start] == fit.iterations
        @test fit.start_converged[fit.best_start] == fit.converged
        @test fit.start_reasons[fit.best_start] === fit.reason
        @test fit.start_residuals[fit.best_start] == fit.residual
        @test length(fit.start_objectives) == fit.nstart == 3
        @test sum(fit.omega; dims=2) ≈ ones(6,1)
        @test sum(fit.theta; dims=1) ≈ ones(1,2)
        @test issorted(vec(sum(fit.omega; dims=1)); rev=true)
        @test fit.membership == [argmax(row) for row in eachrow(fit.omega)]
        @test fit.cluster_sizes == [count(==(i), fit.membership) for i in 1:2]
        @test fit.membership[1] == fit.membership[2] != fit.membership[3] == fit.membership[4]
        @test all(abs.(fit.omega[5:6,:] .- 0.5) .< 0.02)
        _, _, likelihood, objective = gom_test_em(x, fit.theta, fit.omega, fit.shape)
        @test fit.loglikelihood ≈ likelihood atol=1e-10
        @test fit.objective ≈ objective atol=1e-10
        direct = indicators(fit.theta; species_labels=species, region_labels=fit.region_labels)
        linked = indicators(fit)
        @test isequal(linked.top_features, direct.top_features)
        @test isequal(linked.scores, direct.scores)
        @test linked.species_labels == species && linked.region_labels == fit.region_labels
        @test data.matrix == before
    end
    one_start = fitgom(data; k=2, nstart=1, rng=StableRNG(23))
    many_starts = fitgom(data; k=2, nstart=4, rng=StableRNG(23))
    @test many_starts.objective >= one_start.objective
    @test many_starts.start_objectives[1] == one_start.objective
    short = fitgom(data; k=2, nstart=1, maxiter=1, tol=1e-14, rng=StableRNG(42))
    @test !short.converged && short.reason === :maxiter && short.iterations == 1
    # Changes can initially grow away from symmetric factors. A small
    # preliminary update must not imply convergence of the returned factors.
    symmetric_counts = [100 0; 0 100]
    symmetric_initial = [0.5001 0.4999; 0.4999 0.5001]
    premature = fitgom(symmetric_counts; k=2, nstart=1, maxiter=1, tol=3e-4,
        accelerate=false, init_theta=symmetric_initial, init_omega=symmetric_initial)
    @test !premature.converged
    @test premature.reason === :maxiter
    @test premature.residual > premature.tol
    expected_t, expected_w, _, _ = gom_test_em(symmetric_counts,
        premature.theta, premature.omega, premature.shape)
    @test premature.residual ≈ max(maximum(abs.(expected_t .- premature.theta)),
        maximum(abs.(expected_w .- premature.omega))) atol=1e-14
    continued = fitgom(symmetric_counts; k=2, nstart=1, tol=3e-4,
        accelerate=false, init_theta=symmetric_initial, init_omega=symmetric_initial)
    @test continued.converged && continued.reason === :tolerance
    @test continued.residual <= continued.tol
    @test continued.iterations > premature.iterations
    initial_t = [0.4 0.1; 0.3 0.2; 0.2 0.3; 0.1 0.4]
    initial_w = repeat([0.6 0.4], 6, 1)
    rng = StableRNG(9)
    control = copy(rng)
    warm = fitgom(data; k=2, nstart=1, rng, init_theta=initial_t, init_omega=initial_w)
    @test rand(rng) == rand(control)
    @test initial_t == [0.4 0.1; 0.3 0.2; 0.2 0.3; 0.1 0.4]
    @test initial_w == repeat([0.6 0.4], 6, 1)
    warm.theta[1,1] = 0
    warm.sample_labels[1] = "changed"
    @test data.sample_labels == samples && initial_t[1,1] == 0.4
    tight = fitgom(data; k=2, nstart=1, tol=1e-16, accelerate=false, rng=StableRNG(42))
    @test tight.converged == (tight.residual <= tight.tol)
    @test tight.reason in (:stalled, :tolerance, :maxiter)
    three = [20 10 0 0 0 0; 10 20 0 0 0 0; 0 0 20 10 0 0;
             0 0 10 20 0 0; 0 0 0 0 20 10; 0 0 0 0 10 20]
    three_fit = fitgom(three; k=3, nstart=5, rng=StableRNG(206))
    @test three_fit.k == 3 && three_fit.shape == 1/18
    @test three_fit.membership[1] == three_fit.membership[2]
    @test three_fit.membership[3] == three_fit.membership[4]
    @test three_fit.membership[5] == three_fit.membership[6]
    @test length(unique(three_fit.membership)) == 3
end

@testset "GoM: single region, zero columns and aggregation" begin
    for singleton_counts in ([1 0; 0 1], [50 0; 0 50])
        singleton_fit = fitgom(singleton_counts; k=2, nstart=1, rng=StableRNG(12))
        @test all(>(0), singleton_fit.theta) && all(>(0), singleton_fit.omega)
        @test all(<(1), singleton_fit.omega) && isfinite(singleton_fit.objective)
    end
    x = [2 3 0; 6 4 0]
    fit = fitgom(x; k=1, nstart=1, rng=StableRNG(3))
    @test fit.omega == ones(2,1) && fit.membership == [1,1]
    @test fit.theta[:,1] ≈ ([8,7,0] .+ 1/3) ./ 16
    @test fit.converged && fit.cluster_sizes == [2]
    @test fit.species_labels == ["species_1", "species_2", "species_3"]
    @test_throws ArgumentError indicators(fit)
    padded = sparse(x)
    padded[1,3] = 1
    padded[1,3] = 0
    @test nnz(padded) > count(>(0), padded)
    actual = fitgom(padded; k=1, nstart=1, rng=StableRNG(3))
    @test actual.theta == fit.theta
    @test actual.objective == fit.objective
    @test nnz(padded) > count(>(0), padded)
    fraction = fitgom(x ./ 10; k=1, nstart=1, rng=StableRNG(3))
    @test fraction.theta[:,1] ≈ ([0.8,0.7,0] .+ 1/3) ./ 2.5
    @test fraction.theta != fit.theta
    data = CommunityData([1 0 1 0; 0 1 0 1; 1 1 0 0]; species_labels=["A","B","C","D"])
    grouped = aggregate_species(data, ["west", "west", "east", "east"])
    aggregated = fitgom(grouped; k=2, nstart=1, rng=StableRNG(21))
    @test aggregated.sample_labels == data.sample_labels
    @test aggregated.species_labels == grouped.species_labels
    tree = mktempdir() do directory
        path = joinpath(directory, "tree.nwk")
        write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
        read_tree(path)
    end
    lineages = collapse_range(data, tree; k=2).community
    @test fitgom(lineages; k=2, nstart=1, rng=StableRNG(21)).species_labels == lineages.species_labels
end

@testset "GoM: rejected inputs and sparse storage" begin
    x = [1 2; 3 4]
    for bad in (zeros(0,2), zeros(2,0), zeros(2,2), [1 0; 0 0],
        [-1 2; 3 4], [NaN 2; 3 4], [Inf 2; 3 4], [missing 2; 3 4])
        @test_throws ArgumentError fitgom(bad; k=1)
        @test_throws ArgumentError fitgom(sparse(bad); k=1)
    end
    for k in (0, -1, 3, true, 1.5, big(typemax(Int))+1)
        @test_throws ArgumentError fitgom(x; k)
    end
    for bad in (0, -1, true, 1.5, big(typemax(Int))+1)
        @test_throws ArgumentError fitgom(x; k=2, nstart=bad)
        @test_throws ArgumentError fitgom(x; k=2, maxiter=bad)
    end
    for bad in (0, -1, Inf, NaN, true, big"1e1000", big"1e-1000", ones(2,2))
        @test_throws ArgumentError fitgom(x; k=2, shape=bad)
    end
    for bad in (0, -1, 1, Inf, NaN, true, big"1e-1000")
        @test_throws ArgumentError fitgom(x; k=2, tol=bad)
    end
    @test_throws ArgumentError fitgom(x; k=2, accelerate=1)
    @test_throws ArgumentError fitgom(fill(big"1e1000",2,2); k=2)
    @test_throws ArgumentError fitgom(fill(big"1e-1000",2,2); k=2)
    @test_throws ArgumentError fitgom(fill(1e308,2,2); k=2)
    @test_throws ArgumentError fitgom(x; k=2, sample_labels=["same","same"])
    @test_throws ArgumentError fitgom(x; k=2, species_labels=["A"])
    @test_throws ArgumentError fitgom(CommunityData(x); k=2, sample_labels=["a","b"])
    @test_throws ArgumentError fitgom(CommunityData(x); k=2, species_labels=["A","B"])
    for bad in (zeros(2,2), [-1 1; 2 3], [NaN 1; 2 3], [Inf 1; 2 3], ones(3,2), [missing 1; 2 3])
        @test_throws ArgumentError fitgom(x; k=2, init_theta=bad)
        @test_throws ArgumentError fitgom(x; k=2, init_omega=bad)
    end
    n, p = 1000, 20000
    large = sparse(repeat(1:n, 2), vcat(1:n, (p-n+1):p), ones(2n), n, p)
    EvoKits._gom_counts(large)
    allocated = @allocated counts = EvoKits._gom_counts(large)
    @test allocated < 500000
    @test length(counts.values) == 2n && counts.p == p
    @test Base.summarysize(counts) < 100000
    fit = fitgom(large; k=2, nstart=1, maxiter=2, rng=StableRNG(42))
    @test size(fit.theta) == (p,2) && size(fit.omega) == (n,2)
    # Includes 21,000 owned text labels as well as the two factor matrices.
    @test Base.summarysize(fit) < 2_000_000 && isfinite(fit.objective)
end

@testset "GoM: saved local maptpx reference" begin
    reference = CSV.read(joinpath(@__DIR__, "fixtures", "gom", "reference.tsv"), DataFrame; delim='\t')
    for name in unique(reference.case)
        rows = reference[reference.case .== name, :]
        n, p, k = first(rows.n), first(rows.p), first(rows.k)
        x = reshape(parse.(Float64, split(first(rows.counts), ',')), n, p)
        theta = reshape(parse.(Float64, split(first(rows.initial_theta), ',')), p, k)
        omega = reshape(parse.(Float64, split(first(rows.initial_omega), ',')), n, k)
        expected_t = reshape(parse.(Float64, split(first(rows.next_theta), ',')), p, k)
        expected_w = reshape(parse.(Float64, split(first(rows.next_omega), ',')), n, k)
        counts = EvoKits._gom_counts(sparse(x))
        actual_t, actual_w = similar(theta), similar(omega)
        EvoKits._gom_em!(actual_t, actual_w, counts, theta, omega, first(rows.shape), zeros(k))
        @test actual_t ≈ expected_t atol=2e-14
        @test actual_w ≈ expected_w atol=2e-14
        likelihood, objective = EvoKits._gom_score(counts, theta, omega, first(rows.shape), zeros(k))
        @test likelihood ≈ first(rows.loglikelihood) atol=1e-11
        @test objective ≈ first(rows.objective) atol=1e-11
    end
end
