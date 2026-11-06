using SparseArrays, Random

# Reuse the independent original-edge/BigFloat oracle from rao.jl.
function test_rao_summary(actual, expected)
    @test actual isa RaoDResult
    @test size(actual.Dkl) == size(actual.H) == (0, 0)
    @test actual.Dkl !== actual.H
    partition_scale = max(isfinite(expected.total) ? abs(expected.total) : 0,
        isfinite(expected.alpha) ? abs(expected.alpha) : 0)
    for field in (:Dkk, :total, :alpha, :beta, :Fst)
        atol = field === :beta ? 3e-12 * partition_scale :
            field === :Fst && expected.total > 0 ? 3e-12 * partition_scale / expected.total : 0.0
        @test isapprox(getproperty(actual, field), getproperty(expected, field);
            atol, rtol = 3e-12, nans = true)
    end
end

@testset "RaoD summary: independent definition, storage and label alignment" begin
    rng = MersenneTwister(20261109)
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        for trial in 1:10
            n, p = rand(rng, 2:8), rand(rng, 3:12)
            write_tree(path, rao_test_newick(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:9, n, p) ./ 4
            matrix[rand(rng, n, p) .< 0.65] .= 0
            iseven(trial) && (matrix[end, :] .= 0)
            columns, sites = randperm(rng, p), randperm(rng, n)
            lengths = copy(tree.edge_length)
            for weighted in (false, true)
                options = (; abundance_weighted = weighted)
                expected = rao_oracle(matrix, tree; weighted)
                default = raoD(matrix, tree; options...)
                explicit = raoD(matrix, tree; pairwise = true, options...)
                for field in propertynames(default)
                    @test isequal(getproperty(default, field), getproperty(explicit, field))
                end
                for (input, tips) in ((matrix, tree.tip_labels), (sparse(matrix), tree.tip_labels),
                    (view(matrix, :, :), tree.tip_labels), (sparse(matrix[:, columns]), tree.tip_labels[columns]))
                    saved = copy(input)
                    actual = raoD(input, tree; pairwise = false, tip_labels = tips, options...)
                    full = raoD(input, tree; tip_labels = tips, options...)
                    test_rao_summary(actual, expected)
                    for field in (:Dkk, :total, :alpha, :beta, :Fst)
                        @test isequal(getproperty(actual, field), getproperty(full, field))
                    end
                    @test isequal(input, saved)
                end
                data = CommunityData(sparse(matrix[sites, columns]); species_labels = tree.tip_labels[columns])
                test_rao_summary(raoD(data, tree; pairwise = false, options...),
                    rao_oracle(matrix[sites, :], tree; weighted))
            end
            @test tree.edge_length == lengths
        end
        write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
        tree = read_tree(path)
        data = CommunityData(sparse([2 1; 0 3]); species_labels = ["C", "A"])
        actual = raoD(data, tree; pairwise = false)
        test_rao_summary(actual, rao_oracle([1 0 2 0; 3 0 0 0], tree))
        again = raoD(data, tree; pairwise = false)
        @test actual.Dkk !== again.Dkk
        actual.Dkk[1] = -1
        @test again.Dkk[1] > 0
        @test actual.Dkl !== again.Dkl && actual.H !== again.H
    end
end

@testset "RaoD summary: pooling, empty sites and numerical limits" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        for newick in ("(A:1,B:1,C:1);", "(A:0,B:0,C:0);", "A;")
            write_tree(path, newick)
            tree = read_tree(path)
            mixed = zeros(3, tree.ntips)
            mixed[1, :] .= 1
            mixed[2, 1] = 2
            for input in (mixed, sparse(mixed), spzeros(2, tree.ntips), zeros(0, tree.ntips)), weighted in (false, true)
                actual = raoD(input, tree; pairwise = false, abundance_weighted = weighted)
                test_rao_summary(actual, rao_oracle(input, tree; weighted))
            end
        end
        write_tree(path, "((A:1,B:1):1,C:2);")
        tree = read_tree(path)
        cases = (Float16[1 3 0; 0 2 7], Float32[0.1 0.7 0; 0.3 0 0.9],
            BigInt[big(10)^400 big(3)*big(10)^400 0; 0 2 7],
            BigFloat[big"1e400" big"3e400" 0; 0 big"1e-400" big"7e-400"],
            Rational{Int}[1//3 2//3 0; 0 1//7 6//7], Bool[1 1 0; 0 1 1],
            Real[1 3 0; 0 2 7], [1e308 1e308 0; 0 1e308 1e308],
            [1.0 3 0; 0 2 7] .* 1e-300, [1.0 3 0; 0 2 7] .* 1e300,
            [1000 1 0; 1000 1 0], sparse([1, 1, 2], [1, 2, 3], [0.0, 2, 3], 3, 3))
        for matrix in cases, input in (matrix, sparse(matrix)), weighted in (false, true)
            test_rao_summary(raoD(input, tree; pairwise = false, abundance_weighted = weighted),
                rao_oracle(input, tree; weighted))
        end
        incidence = raoD(cases[11], tree; pairwise = false, abundance_weighted = false)
        @test incidence.beta < 0 && incidence.Fst < 0
        for (newick, matrix) in (("((A:1e-100,B:1e-100):1e200,C:1e200);", [2.0 1 0; 0 4 0; 1 0 0]),
            ("(((A:0.5,B:0.5):0.5):20,C:21);", [1.0 1 0; 0 1 2]))
            write_tree(path, newick)
            tree = read_tree(path)
            for weighted in (false, true)
                actual = raoD(sparse(matrix), tree; pairwise = false, abundance_weighted = weighted)
                test_rao_summary(actual, rao_oracle(matrix, tree; weighted))
                @test actual.Dkk[1] > 0 && actual.total > 0
            end
        end
        write_tree(path, "(A:1e307,B:1e307);")
        tree = read_tree(path)
        matrix = [1.0 nextfloat(0.0)]
        actual = raoD(matrix, tree; pairwise = false)
        test_rao_summary(actual, rao_oracle(matrix, tree))
        @test actual.total > 0 && actual.beta == 0
    end
end

@testset "RaoD summary: unchanged input validation" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "(A:1,B:1,C:1);")
        tree = read_tree(path)
        for bad in (-1.0, Inf, NaN), input in ([bad 0 0], sparse([bad 0 0]))
            @test_throws ArgumentError raoD(input, tree; pairwise = false)
        end
        @test_throws ArgumentError raoD(ones(1, 2), tree; pairwise = false)
        for labels in (["A", "B"], ["A", "B", "X"], ["A", "A", "C"])
            @test_throws ArgumentError raoD(ones(1, 3), tree; pairwise = false, tip_labels = labels)
        end
        for bad in (-1.0, Inf, NaN)
            invalid = deepcopy(tree)
            invalid.edge_length[1] = bad
            @test_throws ArgumentError raoD(ones(1, 3), invalid; pairwise = false)
        end
        invalid = deepcopy(tree)
        invalid.edge_length .= floatmax(Float64) / 2
        @test_throws ArgumentError raoD(ones(1, 3), invalid; pairwise = false)
        write_tree(path, "(A:1,B:1.000000005);")
        @test isfinite(raoD(ones(1, 2), read_tree(path); pairwise = false).total)
        write_tree(path, "(A:1,B:1.00000002);")
        @test_throws ArgumentError raoD(ones(1, 2), read_tree(path); pairwise = false)
    end
end

@testset "RaoD summary: many-site storage grows linearly" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        p = 64
        write_tree(path, "(" * join(("sp_$i:1" for i in 1:p), ",") * ");")
        tree = read_tree(path)
        for weighted in (false, true)
            allocations = Int[]
            for n in (500, 6_000)
                input = sparse(repeat(1:n, inner = 3), mod1.(1:3n, p), ones(3n), n, p)
                result = raoD(input, tree; pairwise = false, abundance_weighted = weighted)
                counts = vec(sum(input; dims = 1))
                total = 1 - sum(abs2, counts ./ sum(counts))
                test_rao_summary(result, RaoDResult(Dkk = fill(2 / 3, n), total = total,
                    alpha = 2 / 3, beta = total - 2 / 3, Fst = (total - 2 / 3) / total))
                @test Base.summarysize(result) < 8n + 1_000
                push!(allocations, @allocated raoD(input, tree; pairwise = false, abundance_weighted = weighted))
            end
            @test allocations[2] < 12_000_000
            @test allocations[2] < 16 * allocations[1]
        end
    end
end
