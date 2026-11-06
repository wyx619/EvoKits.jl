using SparseArrays, Random

function rao_virtual_newick(p, shape)
    if shape === :balanced
        parts = ["sp_$i" for i in 1:p]
        while length(parts) > 1
            parts = ["(" * parts[i] * ":1," * parts[i + 1] * ":1)" for i in 1:2:length(parts)]
        end
        return only(parts) * ";"
    elseif shape === :comb
        return repeat("(", p - 1) * "sp_1:1,sp_2:1)" *
            join((":1,sp_$i:$(i - 1))" for i in 3:p)) * ";"
    end
    return "(" * join(("sp_$i:1" for i in 1:p), ",") * ");"
end

@testset "RaoD connecting tree: original-edge oracle and input order" begin
    rng = MersenneTwister(20261110)
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        for trial in 1:16
            p, n = rand(rng, 3:16), rand(rng, 2:7)
            write_tree(path, rao_test_newick(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:8, n, p)
            matrix[rand(rng, n, p) .< 0.8] .= 0
            iseven(trial) && (matrix[end, :] .= 0)
            columns, sites = randperm(rng, p), randperm(rng, n)
            input = sparse(matrix[sites, columns])
            rows = EvoKits._community_distance_rows(input, tree, tree.tip_labels[columns], true;
                need_lengths = false)
            saved = deepcopy(rows)
            original_pool = EvoKits._rao_pool(rows, tree)
            mapped, topology, pooled = EvoKits._rao_prepare(rows, tree, original_pool; strategy = :virtual)
            u = length(unique(reduce(vcat, rows.nodes; init = Int32[])))
            @test u == 0 ? topology === tree : topology.nnodes <= 2u - 1
            @test length(mapped.nodes) == n
            @test all(>=(0), mapped.lengths)
            @test pooled.weights == original_pool.weights && isequal(pooled.site_weights, original_pool.site_weights)
            @test mapped.raw === rows.raw && mapped.weights === rows.weights
            for weighted in (false, true)
                expected = rao_oracle(matrix[sites, :], tree; weighted)
                for strategy in (:full, :virtual, :auto)
                    full = EvoKits._rao_result(rows, tree; strategy, abundance_weighted = weighted)
                    summary = EvoKits._rao_result(rows, tree; strategy, abundance_weighted = weighted, pairwise = false)
                    test_rao_result(full, expected)
                    test_rao_summary(summary, expected)
                    for field in (:Dkk, :total, :alpha, :beta, :Fst)
                        @test isequal(getproperty(summary, field), getproperty(full, field))
                    end
                end
            end
            @test isequal(rows, saved)
            @test input == sparse(matrix[sites, columns])
        end
    end
end

@testset "RaoD connecting tree: automatic sparse and dense selection" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        for shape in (:balanced, :comb, :star)
            p, n = 512, 8
            write_tree(path, rao_virtual_newick(p, shape))
            tree = read_tree(path)
            sparse_input = sparse(repeat(1:n, inner = 3), [mod1(37i, p) for i in 1:3n],
                repeat([1.0, 2, 3], n), n, p)
            # Fully occupied data must preserve the existing full-tree route.
            for matrix in (sparse_input, sparse(ones(2, p))), input in (matrix, Matrix(matrix), view(Matrix(matrix), :, :))
                rows = EvoKits._community_distance_rows(input, tree, nothing, true; need_lengths = false)
                pooled = EvoKits._rao_pool(rows, tree)
                mapped, topology, mapped_pool = EvoKits._rao_prepare(rows, tree, pooled)
                if size(matrix, 1) == n
                    u = length(unique(reduce(vcat, rows.nodes)))
                    @test topology !== tree
                    @test topology.nnodes <= 2u - 1
                    @test mapped_pool.nodes != pooled.nodes
                else
                    @test topology === tree && mapped.nodes === rows.nodes
                end
                @test mapped_pool.weights == pooled.weights
                for weighted in (false, true), pairwise in (false, true)
                    full = EvoKits._rao_result(rows, tree; strategy = :full, abundance_weighted = weighted, pairwise)
                    actual = raoD(input, tree; abundance_weighted = weighted, pairwise)
                    if pairwise
                        test_rao_result(actual, full)
                    else
                        test_rao_summary(actual, full)
                    end
                    # A labelled subset is padded and aligned by CommunityData.
                    columns = reverse(findall(!iszero, vec(sum(matrix; dims = 1))))
                    data = CommunityData(input[:, columns]; species_labels = tree.tip_labels[columns])
                    labeled = raoD(data, tree; abundance_weighted = weighted, pairwise)
                    pairwise ? test_rao_result(labeled, full) : test_rao_summary(labeled, full)
                end
            end
        end
    end
end

@testset "RaoD connecting tree: pooling, numeric limits and empty input" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "((A:1,B:1):1,C:2);")
        tree = read_tree(path)
        cases = (Float16[1 3 0; 0 2 7], Float32[0.1 0.7 0; 0.3 0 0.9],
            BigInt[big(10)^400 big(3)*big(10)^400 0; 0 2 7],
            BigFloat[big"1e400" big"3e400" 0; 0 big"1e-400" big"7e-400"],
            Rational{Int}[1//3 2//3 0; 0 1//7 6//7], Bool[1 1 0; 0 1 1],
            Real[1 3 0; 0 2 7], [1e308 1e308 0; 0 1e308 1e308],
            [1.0 3 0; 0 2 7] .* 1e-300, [1.0 3 0; 0 2 7] .* 1e300,
            [1000 1 0; 1000 1 0], sparse([1, 1, 2], [1, 2, 3], [0.0, 2, 3], 3, 3))
        for matrix in cases, weighted in (false, true)
            expected = rao_oracle(matrix, tree; weighted)
            rows = EvoKits._community_distance_rows(sparse(matrix), tree, nothing, true; need_lengths = false)
            test_rao_result(EvoKits._rao_result(rows, tree; strategy = :virtual, abundance_weighted = weighted), expected)
            test_rao_summary(EvoKits._rao_result(rows, tree; strategy = :virtual, abundance_weighted = weighted, pairwise = false), expected)
        end
        for newick in ("A;", "(A:0,B:0,C:0);", "(((A:0.5,B:0.5):0.5):20,C:21);")
            write_tree(path, newick)
            tree = read_tree(path)
            mono = spzeros(2, tree.ntips)
            mono[1, 1] = 1
            for input in (mono, spzeros(3, tree.ntips), zeros(0, tree.ntips)), weighted in (false, true)
                rows = EvoKits._community_distance_rows(input, tree, nothing, true; need_lengths = false)
                expected = rao_oracle(input, tree; weighted)
                for strategy in (:full, :virtual, :auto)
                    test_rao_result(EvoKits._rao_result(rows, tree; strategy, abundance_weighted = weighted), expected)
                    test_rao_summary(EvoKits._rao_result(rows, tree; strategy, abundance_weighted = weighted, pairwise = false), expected)
                end
            end
        end
        # Automatic compression must retain short terminal paths below huge stems.
        tail = join(("sp_$i:1e200" for i in 3:1024), ",")
        write_tree(path, "((sp_1:1e-100,sp_2:1e-100):1e200," * tail * ");")
        tree = read_tree(path)
        input = sparse([1, 1, 2], [1, 2, 1], [2.0, 1, 1], 2, 1024)
        for weighted in (false, true)
            expected = RaoDResult(Dkk = [weighted ? 4e-100 / 9 : 0.5e-100, 0.0],
                Dkl = weighted ? [4e-100 / 9 1e-100 / 3; 1e-100 / 3 0] : [0.5e-100 0.5e-100; 0.5e-100 0],
                H = weighted ? [0 1e-100 / 9; 1e-100 / 9 0] : [0 0.25e-100; 0.25e-100 0],
                total = 0.375e-100, alpha = weighted ? 1e-100 / 3 : 0.375e-100,
                beta = weighted ? 1e-100 / 24 : 0.0, Fst = weighted ? 1 / 9 : 0.0)
            test_rao_result(raoD(input, tree; abundance_weighted = weighted), expected)
            test_rao_summary(raoD(input, tree; abundance_weighted = weighted, pairwise = false), expected)
        end
        # Tiny positive weights still contribute after multiplying long branches.
        write_tree(path, "(" * join(("sp_$i:1e307" for i in 1:16), ",") * ");")
        tree = read_tree(path)
        input = sparse([1, 1], [1, 2], [1.0, nextfloat(0.0)], 1, 16)
        rows = EvoKits._community_distance_rows(input, tree, nothing, true; need_lengths = false)
        for pairwise in (false, true)
            expected = EvoKits._rao_result(rows, tree; strategy = :full, pairwise)
            actual = EvoKits._rao_result(rows, tree; strategy = :virtual, pairwise)
            @test actual.total > 0 && actual.beta == 0
            pairwise ? test_rao_result(actual, expected) : test_rao_summary(actual, expected)
        end
    end
end

@testset "RaoD connecting tree: original validation and large-tree storage" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        p, n = 20_000, 32
        write_tree(path, rao_virtual_newick(p, :star))
        tree = read_tree(path)
        input = sparse(repeat(1:n, inner = 3), collect(1:3n), ones(3n), n, p)
        rows = EvoKits._community_distance_rows(input, tree, nothing, true; need_lengths = false)
        mapped, topology, pooled = EvoKits._rao_prepare(rows, tree, EvoKits._rao_pool(rows, tree))
        @test topology.nnodes == 3n + 1
        @test Base.summarysize((mapped, topology, pooled)) < 100_000
        for weighted in (false, true), pairwise in (false, true)
            actual = raoD(input, tree; abundance_weighted = weighted, pairwise)
            @test actual.Dkk ≈ fill(2 / 3, n)
            @test actual.total ≈ 1 - 1 / (3n)
            @test actual.alpha ≈ 2 / 3
            @test @allocated(raoD(input, tree; abundance_weighted = weighted, pairwise)) < 4_000_000
        end
        # Invalid, unobserved branches and tips still undergo original checks.
        for pairwise in (false, true), bad in (-1.0, Inf, NaN)
            invalid = deepcopy(tree)
            invalid.edge_length[end] = bad
            @test_throws ArgumentError raoD(input, invalid; pairwise)
            broken = copy(input)
            broken[n, p] = bad
            @test_throws ArgumentError raoD(broken, tree; pairwise)
        end
        broken = deepcopy(tree)
        broken.dist_from_root[broken.tip_ids[end]] += 1
        @test_throws ArgumentError raoD(input, broken; pairwise = false)
        @test_throws ArgumentError EvoKits._rao_result(rows, tree; strategy = :unknown)
    end
end
