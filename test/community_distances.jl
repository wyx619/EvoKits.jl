using SparseArrays, Random, LinearAlgebra

# Independent oracle: sum original branches on the symmetric difference of
# two root paths, then enumerate only present tip pairs with BigFloat weights.
# It does not use root depths, LCAs, production traversals or rerooting.
function community_distance_oracle(matrix, tree; weighted = false, nearest = false, exclude = false)
    paths = [Set{Int}() for _ in 1:tree.ntips]
    edge_of = Dict(tree.child_of_edge[e] => e for e in 1:tree.nedges)
    for (tip, path) in zip(tree.tip_ids, paths)
        node = tip
        while node != tree.root
            edge = edge_of[node]
            push!(path, edge)
            node = tree.parent_of_node[node]
        end
    end
    distances = [sum((BigFloat(tree.edge_length[e]) for e in symdiff(a, b)); init = BigFloat(0))
        for a in paths, b in paths]
    n = size(matrix, 1)
    result = fill(NaN, n, n)
    present = [findall(>(0), matrix[i, :]) for i in 1:n]
    weights = [[weighted ? BigFloat(matrix[i, j]) : BigFloat(1) for j in present[i]] for i in 1:n]
    for i in 1:n, j in 1:i
        (isempty(present[i]) || isempty(present[j])) && continue
        if !nearest
            value = sum(weights[i][a] * weights[j][b] * distances[x, y]
                for (a, x) in enumerate(present[i]), (b, y) in enumerate(present[j])) /
                (sum(weights[i]) * sum(weights[j]))
        else
            function directional(from, to)
                values, used = BigFloat[], BigFloat[]
                for (k, tip) in enumerate(present[from])
                    candidates = [distances[tip, other] for other in present[to] if !exclude || tip != other]
                    isempty(candidates) && continue
                    push!(values, minimum(candidates))
                    push!(used, weights[from][k])
                end
                return isempty(values) ? nothing : sum(values .* used) / sum(used)
            end
            forward, reverse = directional(i, j), directional(j, i)
            if forward === nothing && reverse === nothing
                continue
            elseif forward === nothing
                value = reverse
            elseif reverse === nothing
                value = forward
            elseif weighted
                value = (forward * sum(weights[i]) + reverse * sum(weights[j])) /
                    (sum(weights[i]) + sum(weights[j]))
            else
                value = (forward + reverse) / 2
            end
        end
        result[i, j] = result[j, i] = Float64(value)
    end
    return result
end

distance_values_match(x, y) = size(x) == size(y) &&
    all(isapprox.(x, y; atol = 0, rtol = 2e-12, nans = true))

function distance_test_newick(rng, p)
    parts = ["sp_$i:" * string(rand(rng, [0.0, 0.125, 0.5, 1.0, 3.0])) for i in 1:p]
    while length(parts) > 1
        shuffle!(rng, parts)
        n = rand(rng, 2:min(4, length(parts)))
        children = splice!(parts, 1:n)
        push!(parts, "(" * join(children, ",") * "):" * string(rand(rng, [0.0, 0.25, 2.0])))
    end
    return only(parts) * ";"
end

@testset "Comdist branch sums, diagonals, empty sites and abundance scaling" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "((A:1,B:2):3,(C:4,D:0):2);")
        tree = read_tree(path)
        matrix = [2.0 1 0 0; 0 3 1 0; 0 0 0 0; 4 0 0 0; 0 0 0 1]
        for weighted in (false, true)
            expected = community_distance_oracle(matrix, tree; weighted)
            for input in (matrix, sparse(matrix), view(matrix, :, :))
                saved = copy(input)
                actual = comdist(input, tree; abundance_weighted = weighted)
                @test distance_values_match(actual, expected)
                @test isequal(actual, transpose(actual))
                @test all(isnan, actual[3, :])
                @test actual[4, 4] == 0
                @test input == saved
            end
            @test distance_values_match(comdist(matrix .* 1e300, tree; abundance_weighted = weighted), expected)
            @test distance_values_match(comdist(matrix .* 1e-300, tree; abundance_weighted = weighted), expected)
        end
        @test size(comdist(zeros(0, 4), tree)) == (0, 0)
        @test all(isnan, comdist(spzeros(2, 4), tree))
        # Stored zeros are skipped and input storage is retained verbatim.
        explicit = sparse([1, 1, 2], [1, 2, 4], [0.0, 2, 3], 3, 4)
        @test distance_values_match(comdist(explicit, tree), community_distance_oracle(Matrix(explicit), tree))
        @test nnz(explicit) == 3
        large = BigInt[big(10)^400 big(10)^400 0 0; 0 0 1 0]
        @test distance_values_match(comdist(large, tree; abundance_weighted = true),
            community_distance_oracle(large, tree; weighted = true))
        @test distance_values_match(comdist(Rational{Int}[1//3 2//3 0 0], tree; abundance_weighted = true),
            community_distance_oracle(Rational{Int}[1//3 2//3 0 0], tree; weighted = true))
    end
end

@testset "Comdist randomized sparse input, labels and original-branch oracle" begin
    rng = MersenneTwister(20261106)
    mktempdir() do dir
        path = joinpath(dir, "random.nwk")
        for _ in 1:30
            n, p = rand(rng, 2:9), rand(rng, 3:20)
            write_tree(path, distance_test_newick(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:8, n, p) ./ 4
            matrix[rand(rng, n, p) .< 0.65] .= 0
            matrix[end, :] .= 0
            columns, sites = randperm(rng, p), randperm(rng, n)
            for weighted in (false, true)
                expected = community_distance_oracle(matrix, tree; weighted)
                for (input, tips) in ((matrix, tree.tip_labels), (sparse(matrix), tree.tip_labels),
                    (sparse(matrix[:, columns]), tree.tip_labels[columns]))
                    actual = comdist(input, tree; tip_labels = tips, abundance_weighted = weighted)
                    @test distance_values_match(actual, expected)
                end
                data = CommunityData(sparse(matrix[sites, columns]); species_labels = tree.tip_labels[columns])
                @test distance_values_match(comdist(data, tree; abundance_weighted = weighted), expected[sites, sites])
            end
        end
    end
end

@testset "Comdist retains short paths below long shared stems" begin
    mktempdir() do dir
        path = joinpath(dir, "stem.nwk")
        write_tree(path, "(((A:1e-100,B:2e-100):1e200,C:1):1e100,D:1);")
        tree = read_tree(path)
        matrix = [1.0 0 0 0; 0 1 0 0; 1 1 0 0]
        for weighted in (false, true), input in (matrix, sparse(matrix .* 1e300))
            actual = comdist(input, tree; abundance_weighted = weighted)
            @test actual[1, 2] ≈ 3e-100
            @test actual[3, 3] ≈ 1.5e-100
            @test distance_values_match(actual, community_distance_oracle(Matrix(input), tree; weighted))
        end
    end
end

@testset "Comdist validation and storage on a large sparse tree community" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "(A:1,B:1,C:1);")
        tree = read_tree(path)
        for bad in (-1.0, Inf, NaN), input in ([bad 0 0], sparse([bad 0 0]))
            @test_throws ArgumentError comdist(input, tree)
        end
        @test_throws ArgumentError comdist(ones(1, 2), tree)
        for labels in (["A", "B"], ["A", "B", "X"], ["A", "A", "C"])
            @test_throws ArgumentError comdist(ones(1, 3), tree; tip_labels = labels)
        end
        for bad in (-1.0, Inf, NaN)
            invalid = deepcopy(tree)
            invalid.edge_length[1] = bad
            @test_throws ArgumentError comdist(ones(1, 3), invalid)
        end
        p, n = 20_000, 12
        write_tree(path, "(" * join(("sp_$i:1" for i in 1:p), ",") * ");")
        tree = read_tree(path)
        matrix = sparse(repeat(1:n, inner = 3), collect(1:3n), fill(1.0, 3n), n, p)
        expected = fill(2.0, n, n)
        expected[diagind(expected)] .= 4 / 3
        for weighted in (false, true)
            @test comdist(matrix, tree; abundance_weighted = weighted) ≈ expected
            @test @allocated(comdist(matrix, tree; abundance_weighted = weighted)) < 15_000_000
        end
    end
end

@testset "Comdistnt directional means, identities and boundary cases" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "((A:0,B:0):1,C:2,D:3);")
        tree = read_tree(path)
        matrix = [10.0 0 0 0; 0 0 1 3; 1 2 0 0; 1 0 0 0; 0 0 0 0]
        for weighted in (false, true), exclude in (false, true), input in (matrix, sparse(matrix))
            saved = copy(input)
            actual = comdistnt(input, tree; abundance_weighted = weighted, exclude_conspecifics = exclude)
            expected = community_distance_oracle(matrix, tree; weighted, nearest = true, exclude)
            @test distance_values_match(actual, expected)
            @test isequal(actual, transpose(actual))
            @test all(isnan, actual[5, :])
            @test input == saved
            @test isequal(actual[1, 4], exclude ? NaN : 0.0)
            @test isequal(actual[1, 1], exclude ? NaN : 0.0)
            # A and B are distinct species separated only by zero-length edges.
            @test actual[3, 3] == 0
            @test actual[3, 1] == 0
            @test distance_values_match(comdistnt(input .* 1e300, tree;
                abundance_weighted = weighted, exclude_conspecifics = exclude), expected)
            @test distance_values_match(comdistnt(input .* 1e-300, tree;
                abundance_weighted = weighted, exclude_conspecifics = exclude), expected)
        end
        write_tree(path, "(A:1,B:2,C:4);")
        tree = read_tree(path)
        # Pin the existing Julia contract where picante merges differently.
        @test comdistnt([1 0 0; 0 1 1], tree)[1, 2] == 3.5
        @test comdistnt([10 0 0; 0 1 3], tree; abundance_weighted = true)[1, 2] ≈ 24 / 7
        # Excluding the dominant conspecific must normalize the remaining
        # eligible weights before they underflow relative to the whole row.
        rare = [1e300 1e-300 2e-300; 1 0 0]
        @test comdistnt(rare, tree; abundance_weighted = true,
            exclude_conspecifics = true)[1, 2] ≈ 13 / 3
        @test distance_values_match(comdistnt(sparse(rare), tree;
            abundance_weighted = true, exclude_conspecifics = true),
            community_distance_oracle(rare, tree; weighted = true, nearest = true, exclude = true))
        @test size(comdistnt(zeros(0, 3), tree)) == (0, 0)
        @test all(isnan, comdistnt(spzeros(2, 3), tree))
        explicit = sparse([1, 1, 2], [1, 2, 3], [0.0, 2, 3], 3, 3)
        @test distance_values_match(comdistnt(explicit, tree),
            community_distance_oracle(Matrix(explicit), tree; nearest = true))
        @test nnz(explicit) == 3
        huge = BigInt[big(10)^400 0 0; 0 1 3]
        @test distance_values_match(comdistnt(huge, tree; abundance_weighted = true),
            community_distance_oracle(huge, tree; weighted = true, nearest = true))
    end
end

@testset "Comdistnt randomized rerooting and original-branch oracle" begin
    rng = MersenneTwister(20261107)
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        for _ in 1:35
            n, p = rand(rng, 2:9), rand(rng, 3:20)
            write_tree(path, distance_test_newick(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:8, n, p) ./ 4
            matrix[rand(rng, n, p) .< 0.65] .= 0
            matrix[end, :] .= 0
            columns, sites = randperm(rng, p), randperm(rng, n)
            for weighted in (false, true), exclude in (false, true)
                expected = community_distance_oracle(matrix, tree; weighted, nearest = true, exclude)
                options = (; abundance_weighted = weighted, exclude_conspecifics = exclude)
                @test distance_values_match(comdistnt(matrix, tree; options...), expected)
                @test distance_values_match(comdistnt(sparse(matrix[:, columns]), tree;
                    tip_labels = tree.tip_labels[columns], options...), expected)
                data = CommunityData(sparse(matrix[sites, columns]); species_labels = tree.tip_labels[columns])
                @test distance_values_match(comdistnt(data, tree; options...), expected[sites, sites])
            end
        end
    end
end

@testset "Comdistnt shared stems, polytomies, validation and sparse storage" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "(((A:1e-100,B:2e-100):1e200,C:1):1e100,D:1);")
        tree = read_tree(path)
        matrix = [1.0 0 0 0; 0 1 0 0; 1 1 0 0]
        for weighted in (false, true), exclude in (false, true)
            actual = comdistnt(sparse(matrix), tree; abundance_weighted = weighted, exclude_conspecifics = exclude)
            @test distance_values_match(actual, community_distance_oracle(matrix, tree; weighted, nearest = true, exclude))
            @test actual[1, 2] ≈ 3e-100
        end
        write_tree(path, "(A:1,B:1,C:1);")
        tree = read_tree(path)
        for bad in (-1.0, Inf, NaN), input in ([bad 0 0], sparse([bad 0 0]))
            @test_throws ArgumentError comdistnt(input, tree)
        end
        @test_throws ArgumentError comdistnt(ones(1, 2), tree)
        for labels in (["A", "B"], ["A", "B", "X"], ["A", "A", "C"])
            @test_throws ArgumentError comdistnt(ones(1, 3), tree; tip_labels = labels)
        end
        for bad in (-1.0, Inf, NaN)
            invalid = deepcopy(tree)
            invalid.edge_length[1] = bad
            @test_throws ArgumentError comdistnt(ones(1, 3), invalid)
        end
        # A large star tests nearest/second-nearest propagation across ties
        # and enforces storage independent of the site-by-species product.
        p, n = 20_000, 12
        write_tree(path, "(" * join(("sp_$i:1" for i in 1:p), ",") * ");")
        tree = read_tree(path)
        matrix = sparse(repeat(1:n, inner = 3), collect(1:3n), fill(1.0, 3n), n, p)
        for weighted in (false, true), exclude in (false, true)
            expected = fill(2.0, n, n)
            expected[diagind(expected)] .= exclude ? 2.0 : 0.0
            @test comdistnt(matrix, tree; abundance_weighted = weighted, exclude_conspecifics = exclude) ≈ expected
            @test @allocated(comdistnt(matrix, tree;
                abundance_weighted = weighted, exclude_conspecifics = exclude)) < 15_000_000
        end
    end
end

@testset "Community distances across numeric types, combs and zero branches" begin
    mktempdir() do dir
        path = joinpath(dir, "edge_cases.nwk")
        for newick in ("A;", "(A:0,B:0,C:0);",
            "((((A:1,B:2):0,C:3):4,D:5):6,E:7);")
            write_tree(path, newick)
            tree = read_tree(path)
            p = tree.ntips
            matrix = zeros(Int, 4, p)
            matrix[1, 1] = 1
            matrix[2, :] .= 1
            matrix[3, p] = 1
            for input in (Bool.(matrix), sparse(Float32.(matrix)), Matrix{Real}(matrix))
                for weighted in (false, true)
                    @test distance_values_match(comdist(input, tree; abundance_weighted = weighted),
                        community_distance_oracle(matrix, tree; weighted))
                    for exclude in (false, true)
                        @test distance_values_match(comdistnt(input, tree;
                            abundance_weighted = weighted, exclude_conspecifics = exclude),
                            community_distance_oracle(matrix, tree; weighted, nearest = true, exclude))
                    end
                end
            end
        end
        write_tree(path, "(A:1e100,B:1e100,C:0);")
        tree = read_tree(path)
        matrix = [1.0 1e-20 0; 1 0 0]
        @test comdist(matrix, tree; abundance_weighted = true)[1, 2] ≈ 2e80
        @test distance_values_match(comdist(matrix, tree; abundance_weighted = true),
            community_distance_oracle(matrix, tree; weighted = true))
        write_tree(path, "((A:1,B:2):3,C:4);")
        tree = read_tree(path)
        for input in (Float16[1 3 0; 0 2 7], Float32[0.1 0.7 0; 0.3 0 0.9],
            BigFloat[big"1e400" big"3e400" 0; 0 big"1e-400" big"7e-400"])
            for stored in (input, sparse(input))
                @test distance_values_match(comdist(stored, tree; abundance_weighted = true),
                    community_distance_oracle(Matrix(stored), tree; weighted = true))
                for exclude in (false, true)
                    @test distance_values_match(comdistnt(stored, tree;
                        abundance_weighted = true, exclude_conspecifics = exclude),
                        community_distance_oracle(Matrix(stored), tree; weighted = true, nearest = true, exclude))
                end
            end
        end
    end
end
