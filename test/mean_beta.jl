using SparseArrays, Random

function posterior_tree(newick)
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, newick)
        read_tree(path)
    end
end
posterior_match(a, b) = size(a) == size(b) && all(isapprox.(a, b; rtol = 3e-12, atol = 0, nans = true))
function posterior_same(a, b)
    @test posterior_match(a.distances, b.distances)
    @test a.counts == b.counts
    @test a.sample_labels == b.sample_labels
    @test a.pairs == b.pairs
    @test a.replicates == b.replicates
    @test a.component == b.component
    @test a.index_family == b.index_family
end

@testset "Direct posterior beta matches component-specific generator means" begin
    rng = MersenneTwister(20261108)
    trees = posterior_tree.(["((A:1,B:2):4,(C:2,D:3):1);",
        "((D:3,B:1):2,(A:1,C:2):0);", "((C:0,A:0):0,(D:0,B:0):0);"])
    for n in (0, 1, 7, 129)
        matrix = rand(rng, 0:3, n, 4)
        n > 0 && (matrix[1, :] .= 0)
        n > 1 && (matrix[2, :] .= [1, 0, 0, 0])
        labels = ["site_$i" for i in 1:n]
        pairs = n == 0 ? Tuple{Int, Int}[] : [(1, n), (n, 1), (n, n), (1, n)]
        for family in (:sorensen, :jaccard), component in (:total, :turnover, :nestedness), policy in (:skip, :propagate)
            options = (; index_family = family, component, nan_policy = policy)
            data = CommunityData(sparse(matrix); sample_labels = labels, species_labels = ["A", "B", "C", "D"])
            for selected in (nothing, pairs)
                baseline = mean_dist((phylobeta(data, tree; index_family = family, pairs = selected) for tree in trees);
                    component, pairs = selected, nan_policy = policy)
                for block in (1, 17, 256)
                    posterior_same(mean_dist(data, (tree for tree in trees); pairs = selected, block_size = block, options...), baseline)
                end
                actual = mean_dist(view(matrix, :, :), trees; sample_labels = labels,
                    tip_labels = ["A", "B", "C", "D"], pairs = selected, options...)
                posterior_same(actual, baseline)
                # Omitted labels mean the first tree's order for every draw.
                order = [findfirst(==(label), ["A", "B", "C", "D"]) for label in trees[1].tip_labels]
                posterior_same(mean_dist(matrix[:, order], trees; sample_labels = labels,
                    pairs = selected, options...), baseline)
            end
        end
    end
end

@testset "Posterior beta preserves long common branches and original labels" begin
    trees = posterior_tree.(["((A:1e-100,B:2e-100):1e200,C:1e200);",
        "((B:4e-100,A:3e-100):1e200,C:1e200);"])
    matrix = [1 0 0; 0 1 0; 1 1 0; 0 0 0]
    data = CommunityData(sparse(matrix); species_labels = ["A", "B", "C"])
    for family in (:sorensen, :jaccard), component in (:total, :turnover, :nestedness)
        expected = mean_dist((phylobeta(data, tree; index_family = family) for tree in trees); component)
        posterior_same(mean_dist(data, trees; index_family = family, component), expected)
        @test mean_dist(data, trees; index_family = family).distances[1, 2] > 0
        # Trigger blocked sparse multiplication with the same extreme paths.
        bigdata = CommunityData(sparse(repeat(matrix, 33, 1)); species_labels = data.species_labels)
        baseline = mean_dist((phylobeta(bigdata, tree; index_family = family) for tree in trees); component)
        posterior_same(mean_dist(bigdata, trees; index_family = family, component, block_size = 13), baseline)
    end
    subset = CommunityData(sparse([1 0; 1 1]); species_labels = ["B", "A"])
    extra = posterior_tree("((X:1,A:2):3,B:4);")
    posterior_same(mean_dist(subset, [trees[1], extra]),
        mean_dist(phylobeta(subset, tree) for tree in [trees[1], extra]))
    for values in (BigInt[big(10)^1000 0; 0 1], Real[1//2 1; 0 0], Float16[1 2; 1 0])
        data = CommunityData(values; species_labels = ["A", "B"])
        posterior_same(mean_dist(data, trees), mean_dist(phylobeta(data, tree) for tree in trees))
    end
end

@testset "Posterior beta validates all input and consumes trees once" begin
    tree = posterior_tree("((A:1,B:1):2,C:1);")
    data = CommunityData([1 2 0; 0 0 0]; species_labels = ["A", "B", "C"])
    @test_throws ArgumentError mean_dist(data, CompactTree[])
    @test_throws ArgumentError mean_dist(data, [tree, 1])
    for option in ((; index_family = :other), (; component = :other), (; nan_policy = :other),
        (; block_size = 0), (; block_size = true), (; block_size = 1.2), (; block_size = big(typemax(Int)) + 1))
        @test_throws ArgumentError mean_dist(data, [tree]; option...)
    end
    for labels in (["A", "B", "A"], ["A", "B", "X"], ["A", "B", ""])
        @test_throws ArgumentError mean_dist(data.matrix, [tree]; tip_labels = labels, pairs = [])
    end
    @test_throws ArgumentError mean_dist(data, [tree]; tip_labels = data.species_labels)
    @test_throws ArgumentError mean_dist(data, [tree]; sample_labels = data.sample_labels)
    @test_throws ArgumentError mean_dist(data.matrix, [tree]; sample_labels = ["x", "x"])
    @test_throws ArgumentError mean_dist(data.matrix, [tree]; pairs = [(true, 1)])
    @test_throws ArgumentError mean_dist(ones(2, 2), [tree])
    for bad in (-1.0, Inf, NaN)
        matrix = Float64.(data.matrix)
        matrix[2, 3] = bad
        @test_throws ArgumentError mean_dist(sparse(matrix), [tree]; pairs = [(1, 1)])
        invalid = deepcopy(tree)
        invalid.edge_length[end] = bad
        @test_throws ArgumentError mean_dist(data, [tree, invalid]; pairs = [])
    end
    invalid = deepcopy(tree)
    invalid.edge_length .= 1e308
    @test_throws ArgumentError mean_dist(data, [invalid])
    @test_throws ArgumentError mean_dist(data, [tree, posterior_tree("(A:1,B:1,X:1);")])
    consumed = Ref(0)
    iterator = (begin
        consumed[] += 1
        consumed[] <= 3 || error("tree iterator revisited")
        tree
    end for _ in 1:3)
    result = mean_dist(data, iterator)
    @test consumed[] == result.replicates == 3
    saved = copy(result.distances)
    fill!(data.matrix, 0)
    fill!(tree.edge_length, 0)
    @test isequal(saved, result.distances)
    @test result.distances !== result.counts
end

@testset "Direct posterior beta keeps selected endpoints and file streams bounded" begin
    tree = posterior_tree("(" * join(("sp_$i:1" for i in 1:40_000), ",") * ");")
    n = 6000
    matrix = sparse(repeat(1:n, inner = 3), 1:3n, ones(3n), n, tree.ntips)
    data = CommunityData(matrix; species_labels = tree.tip_labels)
    pairs = [(1, n), (n, 1), (3, 4), (1, 1)]
    result = mean_dist(data, [tree, tree]; pairs)
    @test result.distances == [1, 1, 1, 0]
    @test result.counts == fill(2, 4)
    @test @allocated(mean_dist(data, [tree, tree]; pairs)) < 18_000_000
    plan = EvoKits._mean_beta_prepare(data.matrix, data.sample_labels, pairs)
    @test sum(length, plan.rows) == 12
    @test length(plan.unique_pairs) == 3
    @test plan.workspace.sites == [1, 3, 4, n]
    mktempdir() do dir
        path = joinpath(dir, "posterior.nwk")
        write(path, "((A:1,B:1):1,C:1);\n((C:2,B:1):1,A:1);\n")
        small = CommunityData([1 1 0; 0 1 1]; species_labels = ["A", "B", "C"])
        streamed = read_trees(path) do source
            mean_dist(small, source)
        end
        baseline = read_trees(path) do source
            mean_dist(phylobeta(small, t) for t in source)
        end
        posterior_same(streamed, baseline)
        @test phyloregion(streamed; k = 2).sample_labels == small.sample_labels
    end
end
