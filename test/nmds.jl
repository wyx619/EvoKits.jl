using Test, EvoKits, Random, StableRNGs, LinearAlgebra, Statistics, SparseArrays, DataFrames, CSV
include("nmds_reference.jl")

@testset "NMDS: weighted min-max isotonic reference and a BigFloat Guttman step" begin
    rng = StableRNG(306)
    for case in 1:12, ties in (:weak, :strong)
        observed = case <= 6 ? Float64.(rand(rng, 0:4, 10)) : rand(rng, 10)
        observed[1] = max(1.0, observed[1])
        initial = randn(rng, 5, 2)
        matrix = nmds_test_matrix(observed, 5)
        fit = nmds(matrix; init = initial, ties, nstart = 1, maxiter = 0)
        reference = nmds_reference_stress(observed, fit.points, ties)
        @test fit.distances ≈ reference.distances atol=1e-13
        @test fit.disparities ≈ reference.fitted atol=1e-13
        @test fit.stress ≈ reference.stress atol=1e-14
        @test fit.dissimilarities == observed
        @test fit.iterations == 0 && length(fit.stress_history) == 1
        @test norm(fit.points) / sqrt(5) ≈ 1
        @test vec(mean(fit.points; dims = 1)) ≈ zeros(2) atol=1e-14
        @test abs(dot(fit.points[:, 1], fit.points[:, 2])) < 1e-13
        @test fit.start_stress == [fit.stress] && fit.best_start == 1
        if fit.stress > 1e-10
            step = nmds_reference_step(observed, fit.points, ties)
            expected = nmds_reference_stress(observed, step, ties)
            actual = nmds(matrix; init = fit.points, ties, nstart = 1,
                maxiter = 1, tol = 0, stress_tol = 0)
            @test actual.distances ≈ expected.distances atol=2e-13
            @test actual.disparities ≈ expected.fitted atol=2e-13
            @test actual.stress ≈ expected.stress atol=2e-14
            @test actual.stress <= fit.stress + 1e-14
        end
    end
end

@testset "NMDS: iterative fit, starts, stopping diagnostics and RNG sequence" begin
    rng = StableRNG(4)
    source = rand(rng, 12, 5)
    matrix = [norm(source[i, :] - source[j, :]) for i in 1:12, j in 1:12]
    initial = nmds(matrix; nstart = 1, maxiter = 0)
    fit = nmds(matrix; nstart = 4, maxiter = 2000, rng = StableRNG(7))
    @test fit.stress < initial.stress - 0.01
    @test fit.converged && fit.reason in (:stress_tol, :relative_tol, :numerical_stall)
    @test all(diff(fit.stress_history) .<= 0)
    @test length(fit.stress_history) == fit.iterations + 1
    @test fit.stress == minimum(fit.start_stress)
    @test fit.best_start == argmin(fit.start_stress)
    @test fit.converged == fit.start_converged[fit.best_start]
    @test fit.iterations == fit.start_iterations[fit.best_start]
    @test fit.reason == fit.start_reason[fit.best_start]
    @test fit.dims == 2 && fit.ties === :weak && fit.component === :distance
    @test fit.index_family === nothing
    @test nmds_reference_stress(fit.dissimilarities, fit.points, :weak).stress ≈ fit.stress atol=1e-13
    short = nmds(matrix; nstart = 2, maxiter = 2000, rng = StableRNG(7))
    @test fit.start_stress[1:2] == short.start_stress
    @test fit.stress <= short.stress
    repeat_fit = nmds(matrix; nstart = 4, maxiter = 2000, rng = StableRNG(7))
    @test fit.points == repeat_fit.points
    @test fit.stress_history == repeat_fit.stress_history
    for init in (:pcoa, :random, rand(rng, 12, 2))
        actual_rng, expected_rng = StableRNG(48), StableRNG(48)
        result = nmds(matrix; init, nstart = 3, maxiter = 5, rng = actual_rng, tol = 0, stress_tol = 0)
        for _ in 1:(init === :random ? 3 : 2)
            rand(expected_rng, 12, 2)
        end
        @test rand(actual_rng, UInt64) == rand(expected_rng, UInt64)
        @test length(result.start_stress) == 3
    end
    limited = nmds(matrix; nstart = 1, maxiter = 1, tol = 0, stress_tol = 0)
    @test !limited.converged && limited.reason === :maxiter && limited.iterations == 1
    unmodified = nmds(matrix; nstart = 1, maxiter = 0, tol = 0, stress_tol = 0)
    @test !unmodified.converged && unmodified.reason === :maxiter
    @test unmodified.start_stress == [unmodified.stress]
end

@testset "NMDS: ties, duplicates, extreme initial units and rank invariance" begin
    rng = StableRNG(91)
    equal = ones(7, 7) - Matrix{Float64}(I, 7, 7)
    weak = nmds(equal; nstart = 1)
    strong = nmds(equal; ties = :strong, nstart = 2, maxiter = 1000, rng = StableRNG(5))
    @test weak.stress == 0 && weak.iterations == 0
    @test strong.stress > 0.1
    @test maximum(strong.disparities) == minimum(strong.disparities)
    simplex = nmds(equal; dims = 6, ties = :strong, nstart = 1)
    @test simplex.stress < 1e-14
    pair = nmds([0.0 3; 3 0]; dims = 1, nstart = 1)
    @test pair.stress == 0 && pair.points[:, 1] ≈ [1, -1]
    @test_throws ArgumentError nmds([0.0 3; 3 0]) # dims must be < nsites.
    observed = Float64.([0, 1, 1, 2, 2, 3, 4, 2, 1, 5, 3, 2, 4, 5, 1])
    matrix = nmds_test_matrix(observed, 6)
    initial = randn(rng, 6, 2)
    for ties in (:weak, :strong)
        base = nmds(matrix; init = initial, ties, nstart = 1, maxiter = 50, tol = 0, stress_tol = 0)
        for transformed in (matrix .* 1e-200, matrix .* 1e200, matrix .^ 3)
            fit = nmds(transformed; init = initial, ties, nstart = 1, maxiter = 50, tol = 0, stress_tol = 0)
            @test fit.points == base.points
            @test fit.stress_history == base.stress_history
        end
        for scale in (1e-250, 1e250)
            fit = nmds(matrix; init = initial .* scale, ties, nstart = 1, maxiter = 0)
            reference = nmds(matrix; init = initial, ties, nstart = 1, maxiter = 0)
            @test fit.distances ≈ reference.distances atol=1e-13
            @test fit.stress ≈ reference.stress atol=1e-14
        end
        nearly_coincident = [0.0 0; 1e-300 0; 1 0; 0 1; 1 1; 2 -1]
        fit = nmds(matrix; init = nearly_coincident, ties, nstart = 1, maxiter = 40)
        @test isfinite(fit.stress) && all(isfinite, fit.points)
        @test all(diff(fit.stress_history) .<= 0)
        @test nmds_reference_stress(observed, fit.points, ties).stress ≈ fit.stress atol=1e-13
    end
    # A duplicated site with zero dissimilarity remains a real regression pair.
    points = [0.0 0; 0 0; 1 0; 0 1; 1 1; 2 1]
    duplicates = [norm(points[i, :] - points[j, :]) for i in 1:6, j in 1:6]
    fit = nmds(duplicates; init = initial, ties = :strong, nstart = 1)
    @test length(fit.disparities) == 15
    @test first(fit.dissimilarities) == 0
    @test nmds_reference_stress(fit.dissimilarities, fit.points, :strong).stress ≈ fit.stress atol=1e-13
    # Near ties remain distinct; vegan's absolute sqrt(eps) rule is not adopted.
    close = nmds_test_matrix([1.0, 1 + 1e-10, 2], 3)
    fit = nmds(close; init = [0.0 0; 1 0; 2 1], nstart = 1, maxiter = 0, ties = :strong)
    @test fit.disparities ≈ nmds_reference_stress(fit.dissimilarities, fit.points, :strong).fitted
end

@testset "NMDS: beta, mean distances, regions, context and functional inputs" begin
    matrix = [1 0 1 0; 0 1 0 1; 1 1 0 0; 0 0 1 1; 1 0 0 1; 0 1 1 0]
    labels = ["west", "east", "north", "south", "middle", "other"]
    beta = beta_diss(sparse(matrix); sample_labels = labels)
    tree = mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "((a:1,b:2):3,(c:2,d:1):2);")
        read_tree(path)
    end
    phylo = phylobeta(matrix, tree; sample_labels = labels)
    for input in (beta, phylo), component in (:total, :turnover, :nestedness)
        # A nestedness matrix can be all zero; define that as a rejected input.
        distance = getproperty(input, component)
        if all(iszero, distance)
            @test_throws ArgumentError nmds(input; component)
            continue
        end
        result = nmds(input; component, nstart = 2, maxiter = 50, rng = StableRNG(2))
        plain = nmds(distance; sample_labels = labels, nstart = 2, maxiter = 50, rng = StableRNG(2))
        @test result.points == plain.points
        @test result.stress == plain.stress
        @test result.sample_labels == labels
        @test result.component === component && result.index_family === :sorensen
    end
    average = mean_dist([beta, beta])
    averaged = nmds(average; nstart = 1)
    @test averaged.points == nmds(beta; nstart = 1).points
    @test averaged.sample_labels == labels && averaged.component === :total
    context = prepare_phyloregion(beta)
    from_context = nmds(context; nstart = 1)
    @test from_context.points == nmds(beta; nstart = 1).points
    @test from_context.sample_labels == labels && from_context.component === :total
    regions = phyloregion(context; k = 3)
    fit = nmds(regions; nstart = 1)
    @test fit.points == nmds(regions.region_distances; nstart = 1).points
    @test fit.sample_labels == ["1", "2", "3"] && fit.component === :total
    traits = DataFrame(species = ["a", "b", "c", "d"], size = [0, 1, 10, 11])
    data = CommunityData(matrix; sample_labels = labels, species_labels = traits.species)
    functional = functional_beta(data, traits; k = 2, nstart = 2, rng = StableRNG(1))
    @test nmds(functional; nstart = 1).points == nmds(functional.beta; nstart = 1).points
    @test_throws ArgumentError nmds(beta_diss(matrix; pairs = [(1, 2)]))
    @test_throws ArgumentError nmds(phylobeta(matrix, tree; pairs = [(1, 2)]))
    @test_throws ArgumentError nmds(mean_dist([beta]; pairs = [(1, 2)]))
    @test_throws ArgumentError nmds(beta_diss(vcat(matrix, zeros(Int, 1, 4))))
    @test_throws ArgumentError nmds(phyloregion(context; k = 1))
    @test_throws ArgumentError nmds(beta; sample_labels = labels)
    @test_throws ArgumentError nmds(average; sample_labels = labels)
    @test_throws ArgumentError nmds(regions; sample_labels = labels[1:3])
    @test_throws ArgumentError nmds(context; sample_labels = labels)
    corrupt = deepcopy(regions)
    resize!(corrupt.sample_labels, 1) # Regional labels are membership IDs, independent of site labels.
    @test nmds(corrupt; nstart = 1).sample_labels == ["1", "2", "3"]
    corrupt.region_distances[1, 2] = corrupt.region_distances[2, 1] = NaN
    @test_throws ArgumentError nmds(corrupt)
    corrupt_context = deepcopy(context)
    corrupt_context.distances[1] = NaN
    @test_throws ArgumentError nmds(corrupt_context)
end

@testset "NMDS: input validation, label order and ownership" begin
    rng = StableRNG(602)
    points = randn(rng, 6, 2)
    matrix = [norm(points[i, :] - points[j, :]) for i in 1:6, j in 1:6]
    initial = randn(rng, 6, 2)
    labels = ["one", "two", "three", "four", "five", "six"]
    original, saved_initial = copy(matrix), copy(initial)
    result = nmds(matrix; init = initial, sample_labels = labels, nstart = 1, maxiter = 20)
    permutation = [4, 1, 6, 2, 3, 5]
    permuted = nmds(matrix[permutation, permutation]; init = initial[permutation, :],
        sample_labels = labels[permutation], nstart = 1, maxiter = 20)
    @test permuted.sample_labels == labels[permutation]
    @test nmds_reference_distances(permuted.points) ≈ nmds_reference_distances(result.points[permutation, :]) atol=1e-12
    @test permuted.stress ≈ result.stress atol=1e-13
    @test matrix == original && initial == saved_initial
    labels[1] = "modified"
    @test result.sample_labels[1] == "one"
    result.points[1] = 100
    result.dissimilarities[1] = 100
    result.distances[1] = 100
    result.disparities[1] = 100
    @test matrix == original && initial == saved_initial
    for options in ((dims = 0,), (dims = true,), (dims = 6,), (dims = 2.0,),
        (nstart = 0,), (nstart = false,), (nstart = 1.5,), (maxiter = -1,),
        (maxiter = true,), (maxiter = 1.5,), (ties = :wrong,), (tol = -1,),
        (tol = Inf,), (tol = NaN,), (tol = 1,), (tol = true,), (stress_tol = -1,),
        (stress_tol = 1,), (stress_tol = NaN,), (stress_tol = false,), (init = :wrong,),
        (init = zeros(6, 2),), (init = ones(6, 2),), (init = zeros(5, 2),),
        (init = fill(NaN, 6, 2),), (sample_labels = ["one"],),
        (sample_labels = ["x", "x", "3", "4", "5", "6"],))
        @test_throws ArgumentError nmds(matrix; nstart = 1, maxiter = 0, options...)
    end
    for bad in (zeros(0, 0), zeros(1, 1), zeros(3, 3), zeros(3, 2),
        [1.0 2 3; 2 0 1; 3 1 0], [0.0 1 2; 2 0 3; 2 3 0],
        [0.0 -1 2; -1 0 3; 2 3 0], [0.0 NaN 2; NaN 0 3; 2 3 0],
        [0.0 Inf 2; Inf 0 3; 2 3 0])
        @test_throws ArgumentError nmds(bad; nstart = 1, maxiter = 0)
    end
    @test nmds(matrix; sample_labels = collect(1:6), nstart = 1).sample_labels == string.(1:6)
    @test nmds(sparse(matrix); init = initial, nstart = 1, maxiter = 20).stress ≈
        nmds(matrix; init = initial, nstart = 1, maxiter = 20).stress
    tiny_initial = BigFloat.(initial)
    tiny_initial[1] = BigFloat("1e-1000")
    @test_throws ArgumentError nmds(matrix; init = tiny_initial, nstart = 1)
    @test_throws ArgumentError nmds(BigFloat.(matrix) .* BigFloat("1e1000"); nstart = 1)
    @test_throws ArgumentError nmds(BigFloat.(matrix) .* BigFloat("1e-1000"); nstart = 1)
end

@testset "NMDS: fixed local vegan monotone fits" begin
    fixture = CSV.read(joinpath(@__DIR__, "fixtures", "nmds", "r_reference.tsv"), DataFrame; delim = '\t')
    decode(text) = reduce(vcat, [permutedims(parse.(Float64, split(row, ','))) for row in split(text, ';')])
    for config in eachrow(fixture)
        matrix, initial = decode(config.input), decode(config.points)
        fit = nmds(matrix; init = initial, dims = config.dims,
            ties = Symbol(config.ties), nstart = 1, maxiter = 0)
        @test fit.distances ≈ parse.(Float64, split(config.distances, ',')) atol=1e-12
        @test fit.disparities ≈ parse.(Float64, split(config.disparities, ',')) atol=1e-12
        @test fit.stress ≈ config.stress atol=1e-13
        @test fit.ties === Symbol(config.ties)
        @test size(fit.points) == (config.n, config.dims)
    end
end

@testset "NMDS: reused iteration storage" begin
    rng = StableRNG(32)
    source = rand(rng, 40, 5)
    matrix = [norm(source[i, :] - source[j, :]) for i in 1:40, j in 1:40]
    nmds(matrix; init = :random, nstart = 2, maxiter = 5, rng = StableRNG(4))
    bytes = @allocated result = nmds(matrix; init = :random, nstart = 4,
        maxiter = 300, tol = 0, stress_tol = 0, rng = StableRNG(4))
    @test bytes < 5_000_000 # Four starts must not retain pair-by-iteration history.
    @test Base.summarysize(result) < 100_000
    @test length(result.stress_history) <= 301
end
