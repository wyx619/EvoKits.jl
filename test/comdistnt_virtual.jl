# Independent original-branch oracle and fixtures are in community_distances.jl.
@testset "Comdistnt selected pairs, directional weighting and validation" begin
    mktempdir() do dir
        path = joinpath(dir, "pairs.nwk")
        write_tree(path, "((A:0,B:0):1,C:2,D:3);")
        tree = read_tree(path)
        matrix = [10.0 0 0 0; 0 0 1 3; 1 2 0 0; 1 0 0 0; 0 0 0 0]
        pairs = [(2, 1), (1, 2), (2, 1), (1, 1), (3, 3), (5, 2), (4, 1), (2, 2)]
        for weighted in (false, true), exclude in (false, true), input in (matrix, sparse(matrix), view(matrix, :, :))
            expected = community_distance_oracle(matrix, tree; weighted, nearest = true, exclude)
            options = (; abundance_weighted = weighted, exclude_conspecifics = exclude)
            actual = comdistnt(input, tree; pairs = (p for p in pairs), options...)
            @test actual isa Vector{Float64}
            @test distance_values_match(actual, [expected[i, j] for (i, j) in pairs])
            @test isequal(actual[1], actual[2]) && isequal(actual[1], actual[3])
            @test actual[5] == 0 # Distinct identities separated by zero branches.
            @test isequal(actual[4], exclude ? NaN : 0.0)
            @test isempty(comdistnt(input, tree; pairs = [], options...))
            @test distance_values_match(comdistnt(input, tree; pairs = [[2, 1]], options...), [expected[2, 1]])
            data = CommunityData(input[:, [4, 2, 1]]; species_labels = ["D", "B", "A"])
            padded = copy(matrix)
            padded[:, 3] .= 0
            reference = community_distance_oracle(padded, tree; weighted, nearest = true, exclude)
            @test distance_values_match(comdistnt(data, tree; pairs, options...), [reference[i, j] for (i, j) in pairs])
        end
        for pairs in ([(0, 1)], [(1, 6)], [(true, 1)], [(1.0, 2)], [(1,)], [1], [(1, 2, 3)])
            @test_throws ArgumentError comdistnt(matrix, tree; pairs)
        end
        for pairs in ([], [(1, 1)])
            for value in (-1.0, Inf, NaN), sparse_input in (false, true)
                bad = copy(matrix)
                bad[5, 4] = value
                @test_throws ArgumentError comdistnt(sparse_input ? sparse(bad) : bad, tree; pairs)
            end
            for labels in (["A", "B", "C", "X"], ["A", "A", "C", "D"])
                @test_throws ArgumentError comdistnt(matrix, tree; pairs, tip_labels = labels)
            end
            for bad in (-1.0, Inf, NaN)
                invalid = deepcopy(tree)
                invalid.edge_length[end] = bad
                @test_throws ArgumentError comdistnt(matrix, invalid; pairs)
            end
            invalid = deepcopy(tree)
            invalid.edge_length .= 1e308
            @test_throws ArgumentError comdistnt(matrix, invalid; pairs)
            @test_throws ArgumentError comdistnt(ones(5, 3), tree; pairs)
        end
        @test isempty(comdistnt(zeros(0, 4), tree; pairs = []))
        @test_throws ArgumentError comdistnt(zeros(0, 4), tree; pairs = [(1, 1)])
        @test all(isnan, comdistnt(spzeros(5, 4), tree; pairs))
        # Default self-pairs are zero and can be averaged with mean_dist.
        selected = [(1, 2), (2, 1), (1, 1), (5, 2)]
        values = comdistnt(matrix, tree; pairs = selected)
        labeled = (; distances = values, sample_labels = ["s$i" for i in 1:5], pairs = selected)
        average = mean_dist([labeled, labeled]; pairs = selected)
        @test distance_values_match(average.distances, values)
        @test average.counts == [2, 2, 2, 0]
    end
end

@testset "Comdistnt virtual and full paths against original-branch oracle" begin
    rng = MersenneTwister(20261109)
    mktempdir() do dir
        path = joinpath(dir, "random.nwk")
        for trial in 1:36
            p, n = rand(rng, 3:25), rand(rng, 3:10)
            write_tree(path, distance_test_newick(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:11, n, p)
            matrix[rand(rng, n, p) .< 0.8] .= 0
            matrix[end, :] .= 0
            columns, sites = randperm(rng, p), randperm(rng, n)
            input = sparse(matrix[sites, columns])
            saved, branches = copy(input), copy(tree.edge_length)
            selected = [(rand(rng, 1:n), rand(rng, 1:n)) for _ in 1:24]
            append!(selected, [(1, 2), (2, 1), (1, 1)])
            for weighted in (false, true), exclude in (false, true)
                expected = community_distance_oracle(matrix[sites, :], tree; weighted, nearest = true, exclude)
                rows = EvoKits._community_distance_rows(input, tree, tree.tip_labels[columns], weighted; need_lengths = false)
                for strategy in (:full, :virtual, :auto)
                    full = EvoKits._comdistnt_adaptive(rows, tree, nothing, weighted, exclude; strategy)
                    pairs = EvoKits._comdistnt_adaptive(rows, tree, selected, weighted, exclude; strategy)
                    @test distance_values_match(full, expected)
                    @test distance_values_match(pairs, [expected[i, j] for (i, j) in selected])
                    @test isequal(pairs[end - 2], pairs[end - 1])
                end
                options = (; abundance_weighted = weighted, exclude_conspecifics = exclude)
                @test distance_values_match(comdistnt(input, tree; tip_labels = tree.tip_labels[columns], options...), expected)
                @test distance_values_match(comdistnt(input, tree; pairs = selected, tip_labels = tree.tip_labels[columns], options...),
                    [expected[i, j] for (i, j) in selected])
                data = CommunityData(input; species_labels = tree.tip_labels[columns])
                @test distance_values_match(comdistnt(data, tree; pairs = selected, options...),
                    [expected[i, j] for (i, j) in selected])
            end
            @test input == saved && nonzeros(input) == nonzeros(saved)
            @test tree.edge_length == branches
        end
    end
end

@testset "Comdistnt contraction preserves identities and extreme scales" begin
    mktempdir() do dir
        path = joinpath(dir, "extreme.nwk")
        for newick in ("A;", "(A:0,B:0,C:0);", "((((A:1,B:2):0,C:3):4,D:5):6,E:7);",
            "((((A:1e-100,B:2e-100):1e-100,C:2e-100):1e200,D:1):1e100,E:1);",
            "((A:1e-100,B:1e-100):1e308,C:1);")
            write_tree(path, newick)
            tree = read_tree(path)
            matrix = zeros(Int, 5, tree.ntips)
            matrix[1, 1] = 1
            matrix[2, 1:min(3, tree.ntips)] .= 1:min(3, tree.ntips)
            matrix[3, min(2, tree.ntips)] = 1
            matrix[4, tree.ntips] = 1
            selected = [(1, 3), (3, 1), (2, 2), (4, 1), (5, 2)]
            for weighted in (false, true), exclude in (false, true), input in (matrix, sparse(Float32.(matrix)), Matrix{Real}(matrix))
                expected = community_distance_oracle(matrix, tree; weighted, nearest = true, exclude)
                rows = EvoKits._community_distance_rows(input, tree, nothing, weighted; need_lengths = false)
                for strategy in (:virtual, :full)
                    @test distance_values_match(EvoKits._comdistnt_adaptive(rows, tree, selected, weighted, exclude; strategy),
                        [expected[i, j] for (i, j) in selected])
                end
                @test distance_values_match(comdistnt(input, tree; pairs = selected,
                    abundance_weighted = weighted, exclude_conspecifics = exclude), [expected[i, j] for (i, j) in selected])
            end
        end
        write_tree(path, "((A:1,B:2):3,C:4);")
        tree = read_tree(path)
        for input in (Float16[1 3 0; 0 2 7], Float32[0.1 0.7 0; 0.3 0 0.9],
            BigInt[big(10)^400 2big(10)^400 0; 0 1 3], Rational{Int}[1//3 2//3 0; 0 1//4 3//4],
            BigFloat[big"1e400" big"3e400" 0; 0 big"1e-400" big"7e-400"], [1e300 1e-300 2e-300; 1 0 0]),
            stored in (input, sparse(input)), exclude in (false, true)
            expected = community_distance_oracle(input, tree; weighted = true, nearest = true, exclude)
            selected = [(2, 1), (1, 2), (1, 1), (2, 2)]
            rows = EvoKits._community_distance_rows(stored, tree, nothing, true; need_lengths = false)
            @test distance_values_match(EvoKits._comdistnt_adaptive(rows, tree, selected, true, exclude; strategy = :virtual),
                [expected[i, j] for (i, j) in selected])
            @test distance_values_match(comdistnt(stored, tree; pairs = selected, abundance_weighted = true,
                exclude_conspecifics = exclude), [expected[i, j] for (i, j) in selected])
        end
    end
end

@testset "Comdistnt selected scheduling and linear sparse storage" begin
    mktempdir() do dir
        path = joinpath(dir, "large.nwk")
        p, n = 40_000, 6_000
        write_tree(path, "(" * join(("sp_$i:1" for i in 1:p), ",") * ");")
        tree = read_tree(path)
        matrix = sparse(repeat(1:n, inner = 3), 1:3n, repeat([1, 2, 3], n), n, p)
        selected = [(1, n), (3, 4), (4, 3), (1, 1)]
        used = falses(n)
        used[[1, 3, 4, n]] .= true
        rows = EvoKits._community_distance_rows(matrix, tree, nothing, true; active_sites = used, need_lengths = false)
        plan = EvoKits._comdistnt_pair_plan(rows, selected)
        @test sum(length, rows.nodes) == 12
        @test plan.npairs == 3 && length(plan.jobs) == 5 && count(plan.used) == 4
        @test plan.mapping[2] == plan.mapping[3]
        virtual_rows, topology = EvoKits._comdist_virtual(rows, tree)
        @test topology.nnodes == 13
        @test Base.summarysize((virtual_rows, topology, plan)) < 5_000_000
        for weighted in (false, true), exclude in (false, true)
            options = (; pairs = selected, abundance_weighted = weighted, exclude_conspecifics = exclude)
            @test comdistnt(matrix, tree; options...) ≈ [2, 2, 2, exclude ? 2 : 0]
            @test @allocated(comdistnt(matrix, tree; options...)) < 10_000_000
        end
        p, n = 1024, 8
        tail = join(("sp_$i:1" for i in 4:p), ",")
        write_tree(path, "(((sp_1:1e-100,sp_2:2e-100):1e-100,sp_3:1e-100):1e200," * tail * ");")
        tree = read_tree(path)
        matrix = sparse(repeat(1:n, inner = 2), repeat([1, 2], n), ones(2n), n, p)
        for exclude in (false, true)
            value = exclude ? 3e-100 : 0.0
            @test comdistnt(matrix, tree; exclude_conspecifics = exclude) ≈ fill(value, n, n)
            @test comdistnt(matrix, tree; pairs = [(2, 1), (1, 1)], exclude_conspecifics = exclude) ≈ [value, value]
            # An unrelated dense endpoint with an empty partner is irrelevant.
            mixed = vcat(matrix[1:2, :], sparse(ones(1, p)), spzeros(1, p))
            pairs = [(1, 2), (2, 1), (3, 4)]
            @test distance_values_match(comdistnt(mixed, tree; pairs, exclude_conspecifics = exclude), [value, value, NaN])
        end
    end
end
