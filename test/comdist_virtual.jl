# Uses the independent original-branch oracle in community_distances.jl.
@testset "Comdist selected pairs, canonical symmetry and validation" begin
    mktempdir() do dir
        path = joinpath(dir, "pairs.nwk")
        write_tree(path, "((A:1,B:2):3,C:4,D:0);")
        tree = read_tree(path)
        matrix = [2.0 1 0 0; 0 3 1 0; 0 0 0 0; 0 0 0 1]
        pairs = [(2, 1), (1, 2), (1, 1), (3, 2), (4, 4), (2, 1)]
        for weighted in (false, true), input in (matrix, sparse(matrix), view(matrix, :, :))
            expected = community_distance_oracle(matrix, tree; weighted)
            selected = comdist(input, tree; abundance_weighted = weighted, pairs = (p for p in pairs))
            @test selected isa Vector{Float64}
            @test distance_values_match(selected, [expected[i, j] for (i, j) in pairs])
            @test isequal(selected[1], selected[2]) && isequal(selected[1], selected[6])
            @test selected[5] == 0
            @test isempty(comdist(input, tree; pairs = Tuple{Int, Int}[]))
            @test comdist(input, tree; pairs = [[1, 2]], abundance_weighted = weighted) ≈ [expected[1, 2]]
            data = CommunityData(input[:, [4, 2, 1]]; species_labels = ["D", "B", "A"])
            padded = copy(matrix)
            padded[:, 3] .= 0
            reference = community_distance_oracle(padded, tree; weighted)
            @test distance_values_match(comdist(data, tree; abundance_weighted = weighted, pairs),
                [reference[i, j] for (i, j) in pairs])
        end
        for pairs in ([(0, 1)], [(1, 5)], [(true, 1)], [(1.0, 2)], [(1,)], [1], [(1, 2, 3)])
            @test_throws ArgumentError comdist(matrix, tree; pairs)
        end
        for pairs in (Tuple{Int, Int}[], [(1, 1)])
            bad = copy(matrix)
            bad[4, 4] = -1
            @test_throws ArgumentError comdist(bad, tree; pairs)
            @test_throws ArgumentError comdist(sparse(bad), tree; pairs)
            @test_throws ArgumentError comdist(matrix, tree; pairs, tip_labels = ["A", "B", "C", "X"])
            invalid = deepcopy(tree)
            invalid.edge_length[1] = NaN
            @test_throws ArgumentError comdist(matrix, invalid; pairs)
        end
        @test isempty(comdist(zeros(0, 4), tree; pairs = []))
        @test_throws ArgumentError comdist(zeros(0, 4), tree; pairs = [(1, 1)])
        @test all(isnan, comdist(spzeros(4, 4), tree; pairs))
        offdiagonal = [(1, 2), (2, 1), (3, 2), (1, 2)]
        values = comdist(matrix, tree; pairs = offdiagonal)
        labeled = (; distances = values, sample_labels = ["s$i" for i in 1:4], pairs = offdiagonal)
        average = mean_dist([labeled, labeled]; pairs = offdiagonal)
        @test distance_values_match(average.distances, values)
        @test average.counts == [2, 2, 0, 2]
    end
end

@testset "Comdist virtual and full trees against original-branch pairs" begin
    rng = MersenneTwister(20261108)
    mktempdir() do dir
        path = joinpath(dir, "virtual.nwk")
        for trial in 1:32
            p, n = rand(rng, 3:22), rand(rng, 3:10)
            write_tree(path, distance_test_newick(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:9, n, p)
            matrix[rand(rng, n, p) .< 0.8] .= 0
            matrix[end, :] .= 0
            columns, sites = randperm(rng, p), randperm(rng, n)
            input = sparse(matrix[sites, columns])
            saved = copy(input)
            selected = [(rand(rng, 1:n), rand(rng, 1:n)) for _ in 1:24]
            append!(selected, [(1, 2), (2, 1), (1, 1)])
            for weighted in (false, true)
                expected = community_distance_oracle(matrix[sites, :], tree; weighted)
                rows = EvoKits._community_distance_rows(input, tree, tree.tip_labels[columns], weighted;
                    need_lengths = false)
                for strategy in (:full, :virtual, :auto)
                    full = EvoKits._comdist_adaptive(rows, tree, nothing; strategy)
                    pairs = EvoKits._comdist_adaptive(rows, tree, selected; strategy)
                    @test distance_values_match(full, expected)
                    @test distance_values_match(pairs, [expected[i, j] for (i, j) in selected])
                    @test isequal(pairs[end - 2], pairs[end - 1])
                end
                virtual_rows, virtual = EvoKits._comdist_virtual(rows, tree)
                u = length(unique(vcat(rows.nodes...)))
                @test virtual.nnodes <= max(0, 2u - 1)
                @test all(>=(0), virtual_rows.lengths)
                @test rows.nodes == EvoKits._community_distance_rows(input, tree, tree.tip_labels[columns], weighted).nodes
            end
            @test input == saved && nonzeros(input) == nonzeros(saved)
        end
    end
end

@testset "Comdist contraction preserves extreme scales and numeric inputs" begin
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
            pairs = [(1, 3), (3, 1), (2, 2), (4, 1), (5, 2)]
            for weighted in (false, true), input in (matrix, sparse(Float32.(matrix)), Matrix{Real}(matrix))
                expected = community_distance_oracle(matrix, tree; weighted)
                rows = EvoKits._community_distance_rows(input, tree, nothing, weighted; need_lengths = false)
                for strategy in (:virtual, :full)
                    @test distance_values_match(EvoKits._comdist_adaptive(rows, tree, pairs; strategy),
                        [expected[i, j] for (i, j) in pairs])
                end
            end
        end
        write_tree(path, "((A:1,B:2):3,C:4);")
        tree = read_tree(path)
        pairs = [(2, 1), (1, 2), (1, 1)]
        for input in (Float16[1 3 0; 0 2 7], BigInt[big(10)^400 2big(10)^400 0; 0 1 3],
            Rational{Int}[1//3 2//3 0; 0 1//4 3//4],
            BigFloat[big"1e400" big"3e400" 0; 0 big"1e-400" big"7e-400"])
            expected = community_distance_oracle(input, tree; weighted = true)
            rows = EvoKits._community_distance_rows(input, tree, nothing, true; need_lengths = false)
            @test distance_values_match(EvoKits._comdist_adaptive(rows, tree, pairs; strategy = :virtual),
                [expected[i, j] for (i, j) in pairs])
        end
    end
end

@testset "Comdist active endpoints and selected output stay sparse" begin
    mktempdir() do dir
        path = joinpath(dir, "large.nwk")
        p, n = 40_000, 6_000
        write_tree(path, "(" * join(("sp_$i:1" for i in 1:p), ",") * ");")
        tree = read_tree(path)
        matrix = sparse(repeat(1:n, inner = 3), 1:3n, repeat([1, 2, 3], n), n, p)
        selected = [(1, n), (3, 4), (4, 3), (1, 1)]
        used = falses(n)
        used[[1, 3, 4, n]] .= true
        rows = EvoKits._community_distance_rows(matrix, tree, nothing, true; active_sites = used,
            need_lengths = false)
        @test sum(length, rows.nodes) == 12
        virtual_rows, topology = EvoKits._comdist_virtual(rows, tree)
        @test topology.nnodes == 13
        @test Base.summarysize((virtual_rows, topology)) < 5_000_000
        for weighted in (false, true)
            result = comdist(matrix, tree; pairs = selected, abundance_weighted = weighted)
            @test result ≈ [2, 2, 2, weighted ? 11 / 9 : 4 / 3]
            @test @allocated(comdist(matrix, tree; pairs = selected, abundance_weighted = weighted)) < 10_000_000
        end
        # All observed tips inside one clade: remove a huge shared stem.
        p, n = 1024, 8
        tail = join(("sp_$i:1" for i in 4:p), ",")
        write_tree(path, "(((sp_1:1e-100,sp_2:2e-100):1e-100,sp_3:1e-100):1e200," * tail * ");")
        tree = read_tree(path)
        matrix = sparse(repeat(1:n, inner = 2), repeat([1, 2], n), ones(2n), n, p)
        rows = EvoKits._community_distance_rows(matrix, tree, nothing, false; need_lengths = false)
        @test EvoKits._comdist_use_virtual(rows, tree, n)
        @test comdist(matrix, tree) ≈ fill(1.5e-100, n, n)
        @test EvoKits._comdist_use_virtual(rows, tree, 1)
        @test comdist(matrix, tree; pairs = [(2, 1)]) ≈ [1.5e-100]
        # An unrelated dense endpoint paired only to an empty site must not
        # enlarge the virtual tree needed by the remaining requests.
        mixed = vcat(matrix[1:2, :], sparse(ones(1, p)), spzeros(1, p))
        pairs = [(1, 2), (2, 1), (3, 4)]
        @test distance_values_match(comdist(mixed, tree; pairs), [1.5e-100, 1.5e-100, NaN])
        # Fully occupied input retains the full-tree path.
        dense_rows = EvoKits._community_distance_rows(ones(2, p), tree, nothing, false; need_lengths = false)
        @test !EvoKits._comdist_use_virtual(dense_rows, tree, 2)
    end
end
