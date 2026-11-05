# Lazy input separates the mean accumulator's storage from a quadratic input.
# Its symmetric values and missing observations do not use production indices.
struct MeanStorageDistances <: AbstractMatrix{Float64}
    n::Int
    draw::Int
end
Base.size(matrix::MeanStorageDistances) = (matrix.n, matrix.n)
function Base.getindex(matrix::MeanStorageDistances, i::Int, j::Int)
    lo, hi = minmax(i, j)
    matrix.draw == 2 && mod(lo + hi, 17) == 0 && return NaN
    i == j && return 0.0
    return (mod(lo + hi + matrix.draw, 97) + 1) / 97 * 10.0^(100 * (mod(hi - lo, 5) - 2))
end

@testset "Mean distances keep bounded scratch and reuse result matrices" begin
    for n in (0, 1, 31, 129, 1024)
        state = EvoKits._mean_dist_state(n, nothing)
        # Two full output arrays plus two triangular Float64 work arrays.
        limit = (8 + sizeof(Int)) * n^2 + 8n * (n + 1) + 1024
        @test Base.summarysize(state) <= limit
        sums, counts = state.sums, state.counts
        EvoKits._mean_dist_finish!(state, nothing, n)
        @test state.sums === sums && state.counts === counts
        @test all(isnan, sums)
        @test all(iszero, counts)
    end
    pairs = [(1, 1), (1, 30_000), (30_000, 1)]
    @test Base.summarysize(EvoKits._mean_dist_state(30_000, pairs)) < 4096
end

@testset "Mean distances packed scratch preserves full public results" begin
    n = 1024
    labels = ["storage_$i" for i in 1:n]
    inputs = [(distances = MeanStorageDistances(n, draw), sample_labels = labels) for draw in 1:2]
    mean_dist(inputs)
    allocated = @allocated actual = mean_dist(inputs)
    # Old four-matrix state alone would exceed 33 MB on 64-bit platforms.
    @test allocated < 27_000_000
    @test Base.summarysize(actual) < 18_000_000
    @test isequal(actual.distances, transpose(actual.distances))
    @test actual.counts == transpose(actual.counts)
    @test actual.replicates == 2
    for (i, j) in ((1, 1), (17, 17), (1024, 1024), (1024, 1), (511, 2),
        (129, 128), (512, 256), (1000, 990), (32, 2))
        values = [input.distances[i, j] for input in inputs]
        finite = filter(isfinite, values)
        expected = Float64(sum(BigFloat.(finite)) / length(finite))
        @test actual.distances[i, j] ≈ expected rtol = 3e-15 atol = 0
        @test actual.counts[i, j] == length(finite)
    end
    # A permuted subset leaves output columns entirely unobserved and checks
    # the packed mapping under both missing-value policies.
    order = [129, 17, 1, 128, 2, 31]
    small_inputs = [(distances = view(MeanStorageDistances(129, draw), order, order),
        sample_labels = labels[order]) for draw in 1:2]
    for policy in (:skip, :propagate)
        full = mean_dist(small_inputs; sample_labels = labels[1:129], nan_policy = policy)
        selected = [(i, j) for j in 1:129 for i in 1:129]
        expected, counts = mean_dist_reference(small_inputs, labels[1:129], selected;
            propagate = policy === :propagate)
        @test all(isapprox.(vec(full.distances), expected; nans = true, rtol = 3e-15, atol = 0))
        @test vec(full.counts) == counts
    end
end

@testset "Direct posterior means use the same reduced full storage" begin
    n, p = 1024, 64
    tree = posterior_tree("(" * join(("sp_$i:1" for i in 1:p), ",") * ");")
    data = CommunityData(sparse(1:n, [1 + mod(i - 1, p) for i in 1:n], ones(n), n, p);
        species_labels = tree.tip_labels)
    mean_dist(data, [tree, tree]; block_size = 17)
    allocated = @allocated actual = mean_dist(data, [tree, tree]; block_size = 17)
    @test allocated < 28_000_000
    @test Base.summarysize(actual) < 18_000_000
    @test isequal(actual.distances, transpose(actual.distances))
    @test all(==(2), actual.counts)
    for (i, j) in ((1, 1), (n, n), (n, 1), (65, 1), (129, 128), (900, 260))
        @test actual.distances[i, j] == (mod(i - j, p) == 0 ? 0.0 : 1.0)
    end
    saved = copy(actual.distances)
    fill!(data.matrix, 0)
    fill!(tree.edge_length, 0)
    @test isequal(actual.distances, saved)
end
