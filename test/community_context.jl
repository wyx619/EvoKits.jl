function check_context_result(actual, expected)
    @test typeof(actual) == typeof(expected)
    for field in propertynames(expected)
        left, right = getproperty(actual, field), getproperty(expected, field)
        if left isa AbstractArray{<:Real}
            @test same_values(left, right)
        else
            @test isequal(left, right)
        end
    end
end

@testset "Prepared community: independent oracle and input order" begin
    rng = MersenneTwister(20261102)
    mktempdir() do dir
        path = joinpath(dir, "context.nwk")
        for trial in 1:25
            ntips = rand(rng, 3:25)
            nsites = trial == 25 ? 128 : rand(rng, 2:15)
            write_tree(path, random_community_tree(rng, ntips))
            tree = read_tree(path)
            comm = Float64.(rand(rng, nsites, ntips) .< 0.3)
            comm[end, :] .= 0
            comm[:, end] .= 0
            perm, rowperm = randperm(rng, ntips), randperm(rng, nsites)
            comm = comm[rowperm, perm]
            tips = tree.tip_labels[perm]
            samples = ["site_$i" for i in rowperm]
            oracle = community_reference(comm, tree; tip_labels = tips)
            for input in (comm, sparse(comm))
                context = prepare_community(input, tree; tip_labels = tips, sample_labels = samples)
                @test context isa CommunityContext
                @test pd(context).sample_labels == samples
                @test weighted_endemism(context).species_occupancy == oracle.species_occupancy
                @test phylo_endemism(context).branch_occupancy == oracle.occupancy
                @test same_values(pd(context).pd, oracle.pd)
                @test same_values(pd(context; include_root = false).pd, oracle.pd_no_root)
                @test same_values(weighted_endemism(context).we, oracle.we)
                @test same_values(phylo_endemism(context).pe, oracle.pe)
                @test same_values(phylo_endemism(context; weighted = false).pe, oracle.strict)
                @test same_values(relative_pd(context).rpd, oracle.rpd)
                @test same_values(relative_pd(context).pd_normalized, oracle.normalized)
                @test same_values(relative_pd(context).pd_equal_normalized, oracle.equal_normalized)
                @test pd(context).richness == oracle.richness
                @test sum(weighted_endemism(context).we) ≈ count(>(0), oracle.species_occupancy)
                @test sum(phylo_endemism(context).pe) ≈ sum(tree.edge_length[oracle.occupancy .> 0])
                for family in (:sorensen, :jaccard)
                    full = phylobeta(context; index_family = family)
                    expected = oracle.beta(family)
                    pairs = [(1, nsites), (nsites, 1), (1, 1), (1, nsites)]
                    selected = phylobeta(context; index_family = family, pairs = (p for p in pairs))
                    @test selected.pairs == pairs
                    @test selected.sample_labels == samples
                    for field in (:total, :turnover, :nestedness)
                        @test same_values(getproperty(full, field), getproperty(expected, field))
                        @test same_values(getproperty(selected, field),
                            [getproperty(expected, field)[i, j] for (i, j) in pairs])
                    end
                    # Splitting pair requests retains the complete collection's row indices.
                    first_batch = phylobeta(context; index_family = family, pairs = pairs[1:2])
                    second_batch = phylobeta(context; index_family = family, pairs = pairs[3:4])
                    @test same_values(vcat(first_batch.total, second_batch.total), selected.total)
                end
                for (f, field) in ((unifrac, :jaccard), (phylosor, :sorensen))
                    expected = oracle.beta(field).total
                    f === phylosor && (expected = 1 .- expected)
                    expected[isnan.(expected)] .= 0
                    @test same_values(f(context), expected)
                    @test same_values(f(context), f(input, tree; tip_labels = tips))
                end
                for f in (pd, phylo_endemism, relative_pd, phylobeta)
                    check_context_result(f(context), f(input, tree; tip_labels = tips, sample_labels = samples))
                end
                check_context_result(weighted_endemism(context), weighted_endemism(input; sample_labels = samples))
                check_context_result(pd(context; include_root = false),
                    pd(input, tree; include_root = false, tip_labels = tips, sample_labels = samples))
                check_context_result(phylo_endemism(context; weighted = false),
                    phylo_endemism(input, tree; weighted = false, tip_labels = tips, sample_labels = samples))
                @test isequal(pd_dataframe(context; include_root = false),
                    pd_dataframe(input, tree; include_root = false, tip_labels = tips, sample_labels = samples))
            end
        end
    end
end

@testset "Prepared community: snapshot and result ownership" begin
    mktempdir() do dir
        path = joinpath(dir, "snapshot.nwk")
        for sparse_input in (false, true)
            write_tree(path, "((A:1,B:1):2,C:3);")
            tree = read_tree(path)
            input = [1.0 1 0; 1 0 1; 0 0 1; 0 0 0]
            sparse_input && (input = sparse(input))
            before = copy(input)
            lengths_before = copy(tree.edge_length)
            tips, samples = ["A", "B", "C"], ["AB", "AC", "C", "empty"]
            context = prepare_community(input, tree; tip_labels = tips, sample_labels = samples)
            @test input == before
            @test tree.edge_length == lengths_before
            if sparse_input
                @test nnz(input) == nnz(before)
            end
            pairs = [(1, 2), (4, 1), (1, 2)]
            calculations = (pd, c -> pd(c; include_root = false), weighted_endemism,
                phylo_endemism, c -> phylo_endemism(c; weighted = false), relative_pd,
                phylobeta, c -> phylobeta(c; pairs = pairs))
            expected = map(f -> f(context), calculations)
            distances, similarities = unifrac(context), phylosor(context)
            # Destroy every original input after preparation, including topology.
            sparse_input ? fill!(nonzeros(input), -1) : fill!(input, -1)
            fill!(tree.edge_length, Inf)
            fill!(tree.parent_of_node, 0)
            fill!(tree.tip_ids, 0)
            fill!(tree.tip_labels, "changed")
            fill!(tips, "changed")
            fill!(samples, "changed")
            for (f, result) in zip(calculations, expected)
                check_context_result(f(context), result)
                # Changing any returned label, diagnostic, or value must not poison reuse.
                modified = f(context)
                for field in propertynames(modified)
                    values = getproperty(modified, field)
                    values isa AbstractArray || continue
                    isempty(values) && continue
                    if eltype(values) <: AbstractString
                        fill!(values, "modified")
                    elseif eltype(values) <: Tuple
                        fill!(values, (2, 3))
                    else
                        fill!(values, 0)
                    end
                end
                check_context_result(f(context), result)
            end
            @test pairs == [(1, 2), (4, 1), (1, 2)]
            uni, sor = unifrac(context), phylosor(context)
            fill!(uni, -1)
            fill!(sor, -1)
            @test same_values(unifrac(context), distances)
            @test same_values(phylosor(context), similarities)
        end
    end
end

@testset "Prepared community: empty, zero, extreme and invalid inputs" begin
    mktempdir() do dir
        path = joinpath(dir, "edges.nwk")
        for newick in ("((A:0,B:0):0,C:0);",
            "((A:0.000001,B:0.000002):1000000000000,C:1);",
            "((A:1e308,B:1e-308):1e-308,C:0);",
            "((A:1e308,B:5e-324):5e-324,C:0);")
            write_tree(path, newick)
            tree = read_tree(path)
            for comm in (zeros(0, 3), zeros(128, 3), [1 0 0; 0 1 0], repeat([1 0 0; 0 1 0], 64, 1))
                context = prepare_community(sparse(comm), tree)
                for f in (pd, phylo_endemism, relative_pd, phylobeta)
                    check_context_result(f(context), f(comm, tree))
                end
                check_context_result(pd(context; include_root = false), pd(comm, tree; include_root = false))
                for family in (:sorensen, :jaccard)
                    pairs = isempty(comm) ? Tuple{Int, Int}[] : [(1, 2), (2, 1)]
                    check_context_result(phylobeta(context; index_family = family, pairs = pairs),
                        phylobeta(comm, tree; index_family = family, pairs = pairs))
                end
                @test same_values(unifrac(context), unifrac(comm, tree))
                @test same_values(phylosor(context), phylosor(comm, tree))
            end
        end
        write_tree(path, "((A:1,B:1):2,C:3);")
        tree = read_tree(path)
        stored_zero = sparse([1, 1, 2], [1, 2, 3], [0.0, 1e-20, 1.0], 4, 3)
        context = prepare_community(stored_zero, tree)
        @test weighted_endemism(context).we == [1, 1, 0, 0]
        @test nnz(stored_zero) == 3
        @test_throws ArgumentError prepare_community(ones(1, 2), tree)
        @test_throws ArgumentError prepare_community(ones(1, 3), tree; tip_labels = ["A", "A", "C"])
        @test_throws ArgumentError prepare_community(ones(1, 3), tree; tip_labels = ["A", "B", "missing"])
        @test_throws ArgumentError prepare_community(ones(1, 3), tree; sample_labels = ["a", "b"])
        for bad in (-1.0, Inf, NaN)
            @test_throws ArgumentError prepare_community([bad 0 0], tree)
            @test_throws ArgumentError prepare_community(sparse([bad 0 0]), tree)
            invalid_tree = deepcopy(tree)
            invalid_tree.edge_length[1] = bad
            @test_throws ArgumentError prepare_community(ones(1, 3), invalid_tree)
        end
        invalid_tree = deepcopy(tree)
        fill!(invalid_tree.edge_length, floatmax(Float64))
        @test_throws ArgumentError prepare_community(ones(1, 3), invalid_tree)
        @test_throws ArgumentError phylobeta(context; index_family = :invalid)
        for bad in ([(0, 1)], [(1, 5)], [(1.0, 2)], [(true, 2)], [(1,)], [1])
            @test_throws ArgumentError phylobeta(context; pairs = bad)
        end
        @test isempty(phylobeta(context; pairs = []).total)
    end
end
