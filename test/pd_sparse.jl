using SparseArrays, Random

function pdtest_tree(newick)
    mktempdir() do dir
        path = joinpath(dir, "pd.nwk")
        write(path, newick)
        read_tree(path)
    end
end

# Independent oracle counts original edges by walking each observed tip to
# the root. It never contracts paths or uses the production tree index.
function pdtest_edges(tree, tips, variant)
    counts = zeros(Int, tree.nnodes)
    for tip in tips
        node = tree.tip_ids[tip]
        while node != tree.root
            counts[node] += 1
            node = tree.parent_of_node[node]
        end
    end
    value = big"0.0"
    for edge in 1:tree.nedges
        child = tree.child_of_edge[edge]
        count = counts[child]
        included = variant === :tbl ? count > 0 :
            variant === :mst ? 0 < count < length(tips) :
            variant === :ueh ? count == 1 :
            variant === :sbl ? count > 1 : count > 0 && tree.is_tip[child]
        included && (value += BigFloat(tree.edge_length[edge]))
    end
    return Float64(value)
end

function pdtest_legacy_sample(probabilities, r, rng)
    keys = [p > 0 ? log(rand(rng)) / p : -Inf for p in probabilities]
    return sort(sortperm(keys; rev = true)[1:r])
end

mutable struct PDTestRNG <: AbstractRNG
    values::Vector{Float64}
    position::Int
end
Random.rand(rng::PDTestRNG) = rng.values[mod1(rng.position += 1, length(rng.values))]

@testset "PD variants against original-edge subset sums" begin
    trees = pdtest_tree.(["((A:1,B:2):4,(C:3,D:0):5);",
        "(((A:1):2,B:0,C:3):4,D:5);", "(A:0,B:0,C:0,D:0);",
        "(((A:1e-100,B:2e-100):1e200,C:3):4,D:5);"])
    for tree in trees
        p = tree.ntips
        matrix = [isodd(mask >> (tip - 1)) ? tip : 0 for mask in 0:(2^p - 1), tip in 1:p]
        for variant in (:tbl, :mst, :ueh, :sbl, :tip)
            expected = [pdtest_edges(tree, findall(>(0), row), variant) for row in eachrow(matrix)]
            for input in (matrix, sparse(matrix), view(matrix, :, :))
                saved = copy(input)
                actual = pd_calc(input, tree; variant)
                @test actual.pd ≈ expected rtol = 3e-13
                @test actual.richness == vec(sum(matrix .> 0; dims = 2))
                @test actual.variant === variant
                @test input == saved
            end
            for row in axes(matrix, 1)
                tips = findall(>(0), matrix[row, :])
                @test pd_calc(tree, tips; variant) ≈ expected[row] rtol = 3e-13
                @test pd_calc(tree, tree.tip_labels[tips]; variant) ≈ expected[row] rtol = 3e-13
            end
            order = [4, 1, 3, 2]
            labels = tree.tip_labels[order]
            @test pd_calc(sparse(matrix[:, order]), tree; tip_labels = labels, variant).pd ≈ expected rtol = 3e-13
            data = CommunityData(sparse(matrix[:, order]); species_labels = labels)
            @test pd_calc(data, tree; variant).pd ≈ expected rtol = 3e-13
        end
    end
    tiny = trees[end]
    @test pd_calc(tiny, [1, 2]; variant = :mst) ≈ 3e-100 rtol = 3e-13
    @test pd_calc(tiny, [1, 2]; variant = :ueh) ≈ 3e-100 rtol = 3e-13
    @test pd_calc(trees[2], "A"; variant = :tip) == 1.0
    @test pd_calc(trees[2], "A"; variant = :ueh) == 7.0
end

@testset "PD sparse paths, dense traversal and input validation" begin
    # species_comb is the iterative fixture shared with species_virtual.jl.
    tree = species_comb(512)
    rng = MersenneTwister(20261109)
    for r in (0, 1, 3, 30, 300, 512), variant in (:tbl, :mst, :ueh, :sbl, :tip)
        tips = sort(randperm(rng, tree.ntips)[1:r])
        @test pd_calc(tree, tips; variant) ≈ pdtest_edges(tree, tips, variant)
    end
    tree = pdtest_tree("((A:1,B:2):3,C:4);")
    for invalid in ([1, 1], [1, "A"], [0], [4], [true], [1.0], ["missing"])
        @test_throws ArgumentError pd_calc(tree, invalid)
    end
    for value in (-1.0, Inf, NaN), input in ([value 0 0], sparse([value 0 0]))
        @test_throws ArgumentError pd_calc(input, tree)
        @test_throws ArgumentError pd_bootstrap(input, tree; reps = 1)
    end
    @test_throws ArgumentError pd_calc([1 0], tree)
    @test_throws ArgumentError pd_calc([1 0 0], tree; tip_labels = ["A", "A", "C"])
    @test_throws ArgumentError pd_calc(tree, [1]; variant = :other)
    for invalid in (0, -1, true, big(typemax(Int)) + 1)
        @test_throws ArgumentError pd_bootstrap(tree, 1; reps = invalid)
        @test_throws ArgumentError pd_bootstrap(tree, invalid; reps = 1)
        @test_throws ArgumentError pd_bootstrap([1 0 0], tree; reps = invalid)
    end
    @test_throws ArgumentError pd_bootstrap(tree, 4; reps = 1)
    for weights in ([0, 0, 0], [1, -1, 0], [1, Inf, 0], [1, NaN, 0], [1, 2])
        @test_throws ArgumentError pd_bootstrap(tree, 1; reps = 1, tip_weights = weights)
    end
    for value in (-1.0, Inf, NaN)
        invalid_tree = deepcopy(tree)
        invalid_tree.edge_length[end] = value
        @test_throws ArgumentError pd_calc(invalid_tree, [1])
        @test_throws ArgumentError pd_bootstrap(invalid_tree, 1; reps = 1)
    end
    huge = [floatmax(Float64), floatmax(Float64), 0.0]
    probabilities = EvoKits._pd_probabilities(tree, huge)
    @test probabilities == [0.5, 0.5, 0.0]
    @test all(isfinite, pd_bootstrap(tree, 1; reps = 8, tip_weights = huge, rng = MersenneTwister(1)))
end

@testset "PD partial sampling preserves full-sort subsets and RNG stream" begin
    sampler = EvoKits._PDSamplingWorkspace()
    for p in (1, 2, 3, 8, 31), r in 1:p
        for raw in (ones(p), Float64.(1:p), [isodd(i) ? 1.0 : 0.0 for i in 1:p])
            probabilities = raw / sum(raw)
            for seed in 1:5
                oldrng, newrng = MersenneTwister(seed), MersenneTwister(seed)
                expected = pdtest_legacy_sample(probabilities, r, oldrng)
                @test EvoKits._pd_sample!(sampler, probabilities, r, newrng) == expected
                @test rand(newrng) == rand(oldrng)
            end
            for sequence in ([0.0], [0.5], [prevfloat(1.0)], [0.0, 0.5, 0.5, prevfloat(1.0)])
                oldrng, newrng = PDTestRNG(sequence, 0), PDTestRNG(sequence, 0)
                @test EvoKits._pd_sample!(sampler, probabilities, r, newrng) ==
                    pdtest_legacy_sample(probabilities, r, oldrng)
                @test newrng.position == oldrng.position
            end
        end
    end
end

@testset "PD bootstrap matches legacy draws and retains boundary semantics" begin
    tree = pdtest_tree("(((A:1,B:2):3,C:4):2,D:5);")
    matrix = [0 0 0 0; 1 0 0 0; 1 1 0 0; 1 1 1 0; 1 1 1 1]
    for weights in (ones(4), [1.0, 2.0, 3.0, 4.0], [0.0, 1.0, 0.0, 2.0]), root in (false, true)
        probabilities = weights / sum(weights)
        oldrng, newrng = MersenneTwister(772), MersenneTwister(772)
        expected = zeros(5, 11)
        for rep in 1:11, site in 2:5
            tips = pdtest_legacy_sample(probabilities, site - 1, oldrng)
            expected[site, rep] = !root && length(tips) == 1 ? NaN :
                pdtest_edges(tree, tips, root ? :tbl : :mst)
        end
        actual = pd_bootstrap(sparse(matrix), tree; reps = 11, tip_weights = weights, include_root = root, rng = newrng)
        @test all(isapprox.(actual.replicates, expected; nans = true, rtol = 1e-13))
        @test rand(oldrng) == rand(newrng)
        @test actual.richness == collect(0:4)
        @test isequal(actual.observed, pd(matrix, tree; include_root = root).pd)
        order = [4, 1, 3, 2]
        data = CommunityData(sparse(matrix[:, order]); species_labels = tree.tip_labels[order])
        @test isequal(pd_bootstrap(data, tree; reps = 11, tip_weights = weights,
            include_root = root, rng = MersenneTwister(772)).replicates, actual.replicates)
    end
    for variant in (:tbl, :mst, :ueh, :sbl, :tip), r in (1, 2, 3, 4)
        oldrng, newrng = MersenneTwister(71), MersenneTwister(71)
        expected = [pdtest_edges(tree, pdtest_legacy_sample(fill(0.25, 4), r, oldrng), variant) for _ in 1:13]
        @test pd_bootstrap(tree, r; reps = 13, method = variant, rng = newrng) ≈ expected
        @test rand(oldrng) == rand(newrng)
    end
    emptyrng = MersenneTwister(19)
    @test pd_bootstrap(zeros(0, 4), tree; reps = 3, rng = emptyrng).replicates == zeros(0, 3)
    @test rand(emptyrng) == rand(MersenneTwister(19))
end

@testset "PD weighted sampling follows sequential without-replacement probabilities" begin
    weights = [1.0, 2.0, 3.0, 4.0]
    probabilities = weights / sum(weights)
    sampler, rng = EvoKits._PDSamplingWorkspace(), MersenneTwister(1909)
    frequencies = zeros(Int, 4, 4)
    draws = 40_000
    for _ in 1:draws
        i, j = EvoKits._pd_sample!(sampler, probabilities, 2, rng)
        frequencies[i, j] += 1
    end
    for i in 1:4, j in (i + 1):4
        expected = probabilities[i] * probabilities[j] / (1 - probabilities[i]) +
            probabilities[j] * probabilities[i] / (1 - probabilities[j])
        @test abs(frequencies[i, j] / draws - expected) < 6sqrt(expected * (1 - expected) / draws)
    end
end

@testset "PD sparse storage and reusable bootstrap working state" begin
    tree = species_comb(20_000)
    matrix = sparse(repeat(1:12; inner = 3), [1 + mod(719i, tree.ntips) for i in 1:36], ones(36), 12, tree.ntips)
    pd_calc(matrix, tree)
    pd_bootstrap(matrix, tree; reps = 2, rng = MersenneTwister(1))
    @test @allocated(pd_calc(matrix, tree)) < 15_000_000
    @test @allocated(pd_bootstrap(matrix, tree; reps = 8, rng = MersenneTwister(1))) < 15_000_000
    probabilities = fill(1 / tree.ntips, tree.ntips)
    sampler = EvoKits._PDSamplingWorkspace()
    EvoKits._pd_sample!(sampler, probabilities, 3, MersenneTwister(1))
    @test length(sampler.keys) == length(sampler.heap) == 3
    @test isempty(sampler.excluded)
    EvoKits._pd_sample!(sampler, probabilities, tree.ntips - 3, MersenneTwister(1))
    @test length(sampler.keys) == length(sampler.heap) == 3
    @test length(sampler.excluded) == tree.ntips
end
