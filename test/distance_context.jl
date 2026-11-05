# The independent original-branch BigFloat oracle is in community_distances.jl.
@testset "Prepared distances: original-branch oracle, modes and labeled inputs" begin
    rng = MersenneTwister(20261110)
    mktempdir() do dir
        path = joinpath(dir, "context.nwk")
        for trial in 1:24
            p, n = rand(rng, 3:22), rand(rng, 3:10)
            write_tree(path, distance_test_newick(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:11, n, p) ./ 4
            matrix[rand(rng, n, p) .< 0.75] .= 0
            matrix[end, :] .= 0
            columns, sites = randperm(rng, p), randperm(rng, n)
            labels = ["site_$i" for i in sites]
            selected = [(rand(rng, 1:n), rand(rng, 1:n)) for _ in 1:12]
            append!(selected, [(1, 2), (2, 1), (1, 1)])
            for input in (matrix[sites, columns], sparse(matrix[sites, columns]))
                context = prepare_community_distances(input, tree;
                    tip_labels = tree.tip_labels[columns], sample_labels = labels)
                fixed = prepare_community_distances(input, tree;
                    tip_labels = tree.tip_labels[columns], sample_labels = labels, pairs = (x for x in selected))
                @test context isa CommunityDistanceContext
                @test context.sample_labels == labels && fixed.sample_labels == labels
                data = CommunityData(input; sample_labels = labels, species_labels = tree.tip_labels[columns])
                labeled = prepare_community_distances(data, tree; pairs = selected)
                for weighted in (false, true)
                    expected = community_distance_oracle(matrix[sites, :], tree; weighted)
                    vector = [expected[i, j] for (i, j) in selected]
                    @test distance_values_match(comdist(context; abundance_weighted = weighted), expected)
                    @test distance_values_match(comdist(context; abundance_weighted = weighted, pairs = selected), vector)
                    @test distance_values_match(comdist(fixed; abundance_weighted = weighted), vector)
                    @test distance_values_match(comdist(labeled; abundance_weighted = weighted, pairs = copy(selected)), vector)
                    for exclude in (false, true)
                        options = (; abundance_weighted = weighted, exclude_conspecifics = exclude)
                        expected = community_distance_oracle(matrix[sites, :], tree; weighted, nearest = true, exclude)
                        vector = [expected[i, j] for (i, j) in selected]
                        @test distance_values_match(comdistnt(context; options...), expected)
                        @test distance_values_match(comdistnt(context; pairs = (x for x in selected), options...), vector)
                        @test distance_values_match(comdistnt(fixed; options...), vector)
                        @test distance_values_match(comdistnt(labeled; options...), vector)
                        @test isequal(comdistnt(fixed; options...)[end - 2], comdistnt(fixed; options...)[end - 1])
                    end
                end
                # New requests remain in the original site index space.
                @test isempty(comdist(context; pairs = []))
                @test isempty(comdistnt(context; pairs = []))
            end
        end
    end
end

@testset "Prepared distances: snapshots, output ownership and shared workspaces" begin
    mktempdir() do dir
        path = joinpath(dir, "snapshot.nwk")
        for p in (4, 512), sparse_input in (false, true)
            tail = join(("sp_$i:1" for i in 3:p), ",")
            write_tree(path, "((sp_1:1e-100,sp_2:2e-100):1e200," * tail * ");")
            tree = read_tree(path)
            input = zeros(4, p)
            input[1, 1:2] .= [2, 1]
            input[2, 1:2] .= [1, 3]
            input[3, 2] = 1
            sparse_input && (input = sparse(input))
            saved = copy(input)
            tips, labels = copy(tree.tip_labels), ["a", "b", "singleton", "empty"]
            pairs = [(1, 2), (2, 1), (1, 1), (4, 3), (3, 1)]
            contexts = (prepare_community_distances(input, tree; tip_labels = tips, sample_labels = labels),
                prepare_community_distances(input, tree; tip_labels = tips, sample_labels = labels, pairs))
            for context in contexts
                functions = (c -> comdist(c; abundance_weighted = true, pairs),
                    c -> comdist(c; pairs), c -> comdistnt(c; pairs),
                    c -> comdistnt(c; abundance_weighted = true, exclude_conspecifics = true, pairs),
                    c -> comdistnt(c; exclude_conspecifics = true, pairs))
                expected = map(f -> f(context), functions)
                @test input == saved
                @test contexts[1].rows.nodes === contexts[1].presence.nodes
                @test contexts[1].rows.lengths === contexts[1].presence.lengths
                @test p == 512 ? context.topology.nnodes == 3 : context.topology.nnodes == tree.nnodes
                # Interleaving both metrics/modes cannot poison their workspaces.
                for _ in 1:3, (f, reference) in zip(functions, expected)
                    actual = f(context)
                    @test distance_values_match(actual, reference)
                    fill!(actual, -1)
                    @test distance_values_match(f(context), reference)
                end
                jobs = [Threads.@spawn(f(context)) for f in functions]
                for (job, reference) in zip(jobs, expected)
                    @test distance_values_match(fetch(job), reference)
                end
            end
            references = [(comdist(c; abundance_weighted = true, pairs),
                comdistnt(c; abundance_weighted = true, exclude_conspecifics = true, pairs)) for c in contexts]
            # Destroy original values, labels and topology, including child lists.
            sparse_input ? fill!(nonzeros(input), -1) : fill!(input, -1)
            fill!(tree.edge_length, Inf)
            fill!(tree.parent_of_node, 0)
            fill!(tree.tip_ids, 0)
            fill!(tree.postorder, 0)
            fill!(tree.preorder, 0)
            foreach(empty!, tree.children)
            fill!(tree.tip_labels, "changed")
            fill!(tips, "changed")
            fill!(labels, "changed")
            fill!(pairs, (4, 4))
            for (c, (mean, nearest)) in zip(contexts, references)
                original_pairs = [(1, 2), (2, 1), (1, 1), (4, 3), (3, 1)]
                @test distance_values_match(comdist(c; abundance_weighted = true, pairs = original_pairs), mean)
                @test distance_values_match(comdistnt(c; abundance_weighted = true,
                    exclude_conspecifics = true, pairs = original_pairs), nearest)
                @test c.sample_labels == ["a", "b", "singleton", "empty"]
            end
        end
    end
end

@testset "Prepared distances: numeric types, scaling and absent species" begin
    mktempdir() do dir
        path = joinpath(dir, "types.nwk")
        write_tree(path, "((A:1,B:2):3,C:4);")
        tree = read_tree(path)
        pairs = [(2, 1), (1, 2), (1, 1), (2, 2)]
        for input in (Bool[1 1 0; 1 0 0], Float16[1 3 0; 0 2 7],
            Float32[0.1 0.7 0; 0.3 0 0.9], Matrix{Real}([1 2 0; 0 3 4]),
            Rational{Int}[1//3 2//3 0; 0 1//4 3//4],
            BigInt[big(10)^400 2big(10)^400 0; 0 1 3],
            BigFloat[big"1e400" big"3e400" 0; 0 big"1e-400" big"7e-400"],
            [1e300 1e-300 2e-300; 1 0 0]), stored in (input, sparse(input))
            context = prepare_community_distances(stored, tree; pairs)
            for weighted in (false, true)
                expected = community_distance_oracle(input, tree; weighted)
                @test distance_values_match(comdist(context; abundance_weighted = weighted), [expected[i, j] for (i, j) in pairs])
                for exclude in (false, true)
                    expected = community_distance_oracle(input, tree; weighted, nearest = true, exclude)
                    @test distance_values_match(comdistnt(context; abundance_weighted = weighted, exclude_conspecifics = exclude),
                        [expected[i, j] for (i, j) in pairs])
                end
            end
            if eltype(input) in (BigInt, BigFloat)
                @test context.rows.raw[1][1] !== input[1, 1]
                @test context.rows.scales[1] !== maximum(input[1, :])
            end
        end
        data = CommunityData(sparse([2 1; 0 3]); species_labels = ["B", "A"])
        context = prepare_community_distances(data, tree)
        padded = [1 2 0; 3 0 0]
        for metric in (comdist, comdistnt)
            @test distance_values_match(metric(context; abundance_weighted = true), metric(padded, tree; abundance_weighted = true))
        end
        for newick in ("A;", "(A:0,B:0,C:0);", "((A:1e-100,B:1e-100):1e308,C:1);")
            write_tree(path, newick)
            tree = read_tree(path)
            matrix = zeros(3, tree.ntips)
            matrix[1, 1] = 1
            matrix[2, :] .= 1
            context = prepare_community_distances(matrix, tree)
            for weighted in (false, true)
                @test distance_values_match(comdist(context; abundance_weighted = weighted),
                    community_distance_oracle(matrix, tree; weighted))
                for exclude in (false, true)
                    @test distance_values_match(comdistnt(context; abundance_weighted = weighted, exclude_conspecifics = exclude),
                        community_distance_oracle(matrix, tree; weighted, nearest = true, exclude))
                end
            end
        end
    end
end

@testset "Prepared distances: validation, fixed scope and empty requests" begin
    mktempdir() do dir
        path = joinpath(dir, "validation.nwk")
        write_tree(path, "(A:1,B:1,C:1);")
        tree = read_tree(path)
        matrix = [1.0 2 0; 0 1 1; 0 0 0]
        for pairs in (nothing, [], [(1, 1)])
            for bad in (-1.0, Inf, NaN), stored in (false, true)
                invalid = copy(matrix)
                invalid[3, 3] = bad
                @test_throws ArgumentError prepare_community_distances(stored ? sparse(invalid) : invalid, tree; pairs)
            end
            @test_throws ArgumentError prepare_community_distances(matrix[:, 1:2], tree; pairs)
            for labels in (["A", "B"], ["A", "B", "X"], ["A", "A", "C"])
                @test_throws ArgumentError prepare_community_distances(matrix, tree; pairs, tip_labels = labels)
            end
            for labels in (["a"], ["a", "a", "b"], ["a", "b", ""])
                @test_throws ArgumentError prepare_community_distances(matrix, tree; pairs, sample_labels = labels)
            end
            for bad in (-1.0, Inf, NaN)
                invalid = deepcopy(tree)
                invalid.edge_length[end] = bad
                @test_throws ArgumentError prepare_community_distances(matrix, invalid; pairs)
            end
        end
        invalid_tree = deepcopy(tree)
        invalid_tree.edge_length .= 1e308
        @test_throws ArgumentError prepare_community_distances(matrix, invalid_tree; pairs = [])
        context = prepare_community_distances(matrix, tree)
        for pairs in ([(0, 1)], [(1, 4)], [(true, 1)], [(1.0, 2)], [(1,)], [1], [(1, 2, 3)])
            @test_throws ArgumentError prepare_community_distances(matrix, tree; pairs)
            @test_throws ArgumentError comdist(context; pairs)
            @test_throws ArgumentError comdistnt(context; pairs)
        end
        fixed = prepare_community_distances(matrix, tree; pairs = [(1, 2), (2, 1), (1, 1)])
        for pairs in (nothing, [], [(1, 3)], [(1, 2)])
            @test_throws ArgumentError comdist(fixed; pairs)
            @test_throws ArgumentError comdistnt(fixed; pairs)
        end
        @test distance_values_match(comdist(fixed), comdist(matrix, tree; pairs = [(1, 2), (2, 1), (1, 1)]))
        @test distance_values_match(comdistnt(fixed), comdistnt(matrix, tree; pairs = [(1, 2), (2, 1), (1, 1)]))
        for empty_input in (zeros(0, 3), zeros(3, 3)), pairs in (nothing, [])
            empty_context = prepare_community_distances(sparse(empty_input), tree; pairs)
            for metric in (comdist, comdistnt)
                expected = pairs === nothing ? fill(NaN, size(empty_input, 1), size(empty_input, 1)) : Float64[]
                @test distance_values_match(metric(empty_context), expected)
            end
        end
        empty_context = prepare_community_distances(matrix, tree; pairs = [])
        @test isempty(comdist(empty_context)) && isempty(comdistnt(empty_context))
        @test_throws ArgumentError comdist(empty_context; pairs = [(1, 1)])
        data = CommunityData(matrix)
        @test_throws ArgumentError prepare_community_distances(data, tree; tip_labels = tree.tip_labels)
        @test_throws ArgumentError prepare_community_distances(data, tree; sample_labels = ["a", "b", "c"])
    end
end

@testset "Prepared distances: fixed endpoints and repeated allocation bounds" begin
    mktempdir() do dir
        path = joinpath(dir, "large.nwk")
        p, n = 40_000, 6_000
        write_tree(path, "(" * join(("sp_$i:1" for i in 1:p), ",") * ");")
        tree = read_tree(path)
        matrix = sparse(repeat(1:n, inner = 3), 1:3n, repeat([1.0, 2, 3], n), n, p)
        pairs = [(1, n), (3, 4), (4, 3), (1, 1)]
        context = prepare_community_distances(matrix, tree; pairs)
        @test sum(length, context.rows.nodes) == 12
        @test context.topology.nnodes == 13
        @test length(context.mean_workspace.mass) == 13 && length(context.nearest_workspace.first_tip) == 13
        @test Base.summarysize(context) < 6_000_000
        @test @allocated(prepare_community_distances(matrix, tree; pairs)) < 15_000_000
        for weighted in (false, true)
            @test comdist(context; abundance_weighted = weighted) ≈ [2, 2, 2, weighted ? 11 / 9 : 4 / 3]
            @test @allocated(comdist(context; abundance_weighted = weighted)) < 5_000
            for exclude in (false, true)
                @test comdistnt(context; abundance_weighted = weighted, exclude_conspecifics = exclude) ≈ [2, 2, 2, exclude ? 2 : 0]
                @test @allocated(comdistnt(context; abundance_weighted = weighted, exclude_conspecifics = exclude)) < 5_000
            end
        end
        # A dense endpoint paired only with an empty site retains no records.
        mixed = vcat(matrix[1:2, :], sparse(ones(1, p)), spzeros(1, p))
        context = prepare_community_distances(mixed, tree; pairs = [(1, 2), (3, 4)])
        @test sum(length, context.rows.nodes) == 6
        @test sum(length, context.rows.raw) == 6 && sum(length, context.presence.raw) == 6
        @test context.topology.nnodes == 7
        @test Base.summarysize(context) < 20_000
        @test distance_values_match(comdist(context), [2, NaN])
        @test distance_values_match(comdistnt(context), [2, NaN])
        view_context = prepare_community_distances(view(matrix, [1, 3, 4, n], :), tree;
            pairs = [(1, 2), (2, 1), (4, 4)])
        @test comdist(view_context) ≈ [2, 2, 4 / 3]
        @test comdistnt(view_context) ≈ [2, 2, 0]
    end
end
