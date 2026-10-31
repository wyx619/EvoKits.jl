using SparseArrays
include("community_reference.jl")

same_values(x, y) = size(x) == size(y) &&
    all(isapprox.(x, y; atol = 1e-12, rtol = 1e-10, nans = true))

@testset "Endemism, relative PD and phylogenetic beta: hand calculations" begin
    mktempdir() do dir
        path = joinpath(dir, "hand.nwk")
        write_tree(path, "((A:1,B:1):2,C:3);")
        tree = read_tree(path)
        comm = [1 1 0; 1 0 1; 0 0 1; 0 0 0]
        labels = ["AB", "AC", "C", "empty"]
        for input in (comm, sparse(comm), Float64.(comm), comm .> 0)
            we = weighted_endemism(input; sample_labels = labels)
            pe = phylo_endemism(input, tree; sample_labels = labels)
            strict = phylo_endemism(input, tree; weighted = false)
            rpd = relative_pd(input, tree; sample_labels = labels)
            @test we.we == [1.5, 1, 0.5, 0]
            @test we.species_occupancy == [2, 1, 2]
            @test we.richness == [2, 2, 1, 0]
            @test we.sample_labels == pe.sample_labels == rpd.sample_labels == labels
            @test pe.pe == [2.5, 3, 1.5, 0]
            @test strict.pe == [1, 0, 0, 0]
            @test strict.weighted == false
            @test pe.weighted == true
            @test pd(input, tree).pd == [4, 6, 3, 0]
            @test same_values(pd(input, tree; include_root = false).pd, [2, 6, NaN, 0])
            @test same_values(rpd.rpd, [16 / 21, 8 / 7, 12 / 7, NaN])
            @test rpd.pd_normalized ≈ [4, 6, 3, 0] ./ 7
            @test rpd.pd_equal_normalized ≈ [3, 3, 1, 0] ./ 4
            @test rpd.total_branch_length == 7
            @test rpd.positive_branch_count == 4
            sor = phylobeta(input, tree)
            jac = phylobeta(input, tree; index_family = :jaccard)
            @test sor.total[1, 2] ≈ 0.4
            @test sor.turnover[1, 2] ≈ 0.25
            @test sor.nestedness[1, 2] ≈ 0.15
            @test jac.total[1, 2] ≈ 4 / 7
            @test jac.turnover[1, 2] ≈ 0.4
            @test jac.nestedness[1, 2] ≈ 6 / 35
            @test sor.total[1, 1] == sor.turnover[1, 1] == sor.nestedness[1, 1] == 0
            @test isnan(sor.total[4, 4])
            @test sor.total[1, 4] == 1
            @test isnan(sor.turnover[1, 4]) && isnan(sor.nestedness[1, 4])
            selected = [(1, 2), (3, 1), (4, 4), (1, 2)]
            sub = phylobeta(input, tree; pairs = (p for p in selected))
            @test sub.pairs == selected
            for field in (:total, :turnover, :nestedness)
                @test same_values(getproperty(sub, field), [getproperty(sor, field)[i, j] for (i, j) in selected])
            end
        end
        explicit_zero = sparse([1, 1, 2], [1, 2, 3], [0.0, 1e-20, 1.0], 4, 3)
        before = copy(explicit_zero)
        @test weighted_endemism(explicit_zero).we == [1, 1, 0, 0]
        @test nnz(explicit_zero) == nnz(before)
        @test nonzeros(explicit_zero) == nonzeros(before)
        @test phylo_endemism(explicit_zero, tree).pe == [3, 3, 0, 0]
        @test same_values(phylobeta(explicit_zero, tree).total, phylobeta(Matrix(explicit_zero), tree).total)
        lengths_before, depth_before = copy(tree.edge_length), copy(tree.dist_from_root)
        relative_pd(comm, tree)
        @test tree.edge_length == lengths_before && tree.dist_from_root == depth_before
        @test size(phylobeta(zeros(0, 3), tree).total) == (0, 0)
        @test isempty(phylobeta(comm, tree; pairs = []).total)
        @test isempty(weighted_endemism(zeros(0, 3)).we)
        @test weighted_endemism(zeros(2, 0)).we == [0, 0]
        @test isempty(phylo_endemism(zeros(0, 3), tree).pe)
        @test isempty(relative_pd(zeros(0, 3), tree).rpd)
        zero_tree = deepcopy(tree)
        fill!(zero_tree.edge_length, 0)
        @test all(isnan, relative_pd(comm, zero_tree).rpd)
        @test all(isnan, phylobeta(comm, zero_tree).total)
        @test all(iszero, phylo_endemism(comm, zero_tree).pe)
        @test all(iszero, unifrac(comm, zero_tree))
        @test all(iszero, phylosor(comm, zero_tree))

        for f in (pd, phylo_endemism, relative_pd, phylobeta, unifrac, phylosor)
            @test_throws ArgumentError f(ones(1, 2), tree)
            @test_throws ArgumentError f(comm, tree; tip_labels = ["A", "B", "X"])
            @test_throws ArgumentError f(comm, tree; tip_labels = ["A", "A", "C"])
            for bad in (-1.0, Inf, NaN)
                @test_throws ArgumentError f([bad 0 0], tree)
                @test_throws ArgumentError f(sparse([bad 0 0]), tree)
            end
            bad_tree = deepcopy(tree)
            bad_tree.edge_length[1] = -1
            @test_throws ArgumentError f(comm, bad_tree)
        end
        for bad in (-1.0, Inf, NaN)
            @test_throws ArgumentError weighted_endemism([bad 0 0])
            @test_throws ArgumentError weighted_endemism(sparse([bad 0 0]))
        end
        for f in (weighted_endemism,)
            @test_throws ArgumentError f(comm; sample_labels = ["wrong"])
        end
        for f in (pd, phylo_endemism, relative_pd, phylobeta)
            @test_throws ArgumentError f(comm, tree; sample_labels = ["wrong"])
        end
        @test_throws ArgumentError phylobeta(comm, tree; index_family = :wrong)
        for bad in ([(0, 1)], [(1, 5)], [(1.0, 2)], [(true, 2)], [(1,)], [1])
            @test_throws ArgumentError phylobeta(comm, tree; pairs = bad)
        end
    end
end

@testset "Sparse community kernels: independent randomized oracle" begin
    rng = MersenneTwister(20261031)
    mktempdir() do dir
        path = joinpath(dir, "random.nwk")
        for trial in 1:100
            ntips, nsites = rand(rng, 3:30), rand(rng, 2:12)
            write_tree(path, random_community_tree(rng, ntips))
            tree = read_tree(path)
            comm = Float64.(rand(rng, nsites, ntips) .< rand(rng))
            trial % 3 == 0 && (comm[1, :] .= 0)
            trial % 7 == 0 && (comm[:, end] .= 0)
            ref = community_reference(comm, tree)
            sparse_comm = sparse(comm)
            for input in (comm, sparse_comm)
                @test same_values(pd(input, tree).pd, ref.pd)
                @test same_values(pd(input, tree; include_root = false).pd, ref.pd_no_root)
                @test same_values(weighted_endemism(input).we, ref.we)
                pe = phylo_endemism(input, tree)
                @test same_values(pe.pe, ref.pe)
                @test pe.branch_occupancy == ref.occupancy
                @test same_values(phylo_endemism(input, tree; weighted = false).pe, ref.strict)
                rpd = relative_pd(input, tree)
                @test same_values(rpd.rpd, ref.rpd)
                @test same_values(rpd.pd_normalized, ref.normalized)
                @test same_values(rpd.pd_equal_normalized, ref.equal_normalized)
                @test sum(weighted_endemism(input).we) ≈ count(>(0), ref.species_occupancy)
                @test sum(pe.pe) ≈ sum(tree.edge_length[ref.occupancy .> 0])
                for family in (:sorensen, :jaccard)
                    beta = phylobeta(input, tree; index_family = family)
                    expected = ref.beta(family)
                    for field in (:total, :turnover, :nestedness)
                        values = getproperty(beta, field)
                        @test same_values(values, getproperty(expected, field))
                        @test same_values(values, transpose(values))
                        @test all(x -> isnan(x) || -1e-12 <= x <= 1 + 1e-12, values)
                    end
                    @test same_values(beta.total[isfinite.(beta.nestedness)],
                        (beta.turnover .+ beta.nestedness)[isfinite.(beta.nestedness)])
                    positive = ref.pd .> 0
                    if family == :jaccard
                        @test same_values(beta.total[positive, positive], unifrac(input, tree)[positive, positive])
                    else
                        @test same_values(beta.total[positive, positive], 1 .- phylosor(input, tree)[positive, positive])
                    end
                end
            end
            perm = randperm(rng, ntips)
            rowperm = randperm(rng, nsites)
            reordered = sparse_comm[rowperm, perm]
            tip_labels = tree.tip_labels[perm]
            @test same_values(pd(reordered, tree; tip_labels = tip_labels).pd, ref.pd[rowperm])
            @test same_values(phylo_endemism(reordered, tree; tip_labels = tip_labels).pe, ref.pe[rowperm])
            @test same_values(relative_pd(reordered, tree; tip_labels = tip_labels).rpd, ref.rpd[rowperm])
            @test same_values(phylobeta(reordered, tree; tip_labels = tip_labels).total,
                ref.beta(:sorensen).total[rowperm, rowperm])
            @test same_values(unifrac(reordered, tree; tip_labels = tip_labels), unifrac(comm, tree)[rowperm, rowperm])
            @test same_values(phylosor(reordered, tree; tip_labels = tip_labels), phylosor(comm, tree)[rowperm, rowperm])
            scaled = deepcopy(tree)
            factor = isodd(trial) ? 1e-12 : 1e12
            scaled.edge_length .*= factor
            scaled.dist_from_root .*= factor
            @test same_values(phylo_endemism(comm, scaled).pe, ref.pe .* factor)
            @test same_values(relative_pd(comm, scaled).rpd, ref.rpd)
            @test same_values(phylobeta(comm, scaled).total, ref.beta(:sorensen).total)
        end
    end
end

@testset "Small unique branches survive large shared branches" begin
    mktempdir() do dir
        path = joinpath(dir, "scale.nwk")
        write_tree(path, "((A:0.000001,B:0.000002):1000000000000,C:1);")
        tree = read_tree(path)
        comm = [1 0 0; 0 1 0]
        beta = phylobeta(comm, tree)
        @test beta.total[1, 2] > 0
        @test beta.turnover[1, 2] > 0
        @test beta.nestedness[1, 2] > 0
        @test beta.total[1, 2] ≈ 1.5e-18
        @test unifrac(comm, tree)[1, 2] ≈ 3e-18
        # Also exercise the full-output sparse multiplication path.
        large = repeat(comm, 64, 1)
        full = phylobeta(sparse(large), tree)
        @test full.total[1, 2] ≈ beta.total[1, 2]
        @test full.turnover[1, 2] ≈ beta.turnover[1, 2]
        @test full.nestedness[1, 2] ≈ beta.nestedness[1, 2]
    end
end

@testset "Beta decomposition across extreme branch scales" begin
    mktempdir() do dir
        path = joinpath(dir, "extreme.nwk")
        for (shared, left, right) in ((1e-308, 1e308, 1e-308),
            (nextfloat(0.0), 1e308, nextfloat(0.0)), (8e307, 1e307, 1.5e307))
            write_tree(path, "((A:$left,B:$right):$shared,C:0);")
            tree = read_tree(path)
            comm = [1 0 0; 0 1 0]
            pairs = [(1, 2), (2, 1), (1, 1), (2, 2)]
            for input in (comm, sparse(comm), sparse(repeat(comm, 64, 1)))
                for family in (:sorensen, :jaccard)
                    # Evaluate the defining fractions with a wider exponent range,
                    # independently of the production scaling and branch traversal.
                    expected = setprecision(BigFloat, 256) do
                        a, b, c = BigFloat.((shared, left, right))
                        m = min(b, c)
                        total = (b + c) / (family == :sorensen ? 2a + b + c : a + b + c)
                        turnover = family == :sorensen ? m / (a + m) : 2m / (a + 2m)
                        Float64.((total, turnover, total - turnover))
                    end
                    full = phylobeta(input, tree; index_family = family, tip_labels = ["A", "B", "C"])
                    selected = phylobeta(input, tree; index_family = family,
                        tip_labels = ["A", "B", "C"], pairs = pairs)
                    for (field, value) in zip((:total, :turnover, :nestedness), expected)
                        @test isapprox(getproperty(full, field)[1, 2], value; atol = 0, rtol = 1e-12)
                        @test isapprox(getproperty(full, field)[2, 1], value; atol = 0, rtol = 1e-12)
                        @test isapprox(getproperty(selected, field), [value, value, 0, 0]; atol = 0, rtol = 1e-12)
                    end
                end
            end
        end
    end
end

@testset "Full-output sparse multiplication matches direct pair calculation" begin
    rng = MersenneTwister(20261101)
    mktempdir() do dir
        path = joinpath(dir, "full.nwk")
        for trial in 1:8
            write_tree(path, random_community_tree(rng, 25))
            tree = read_tree(path)
            comm = sparse(Int.(rand(rng, 128, 25) .< 0.1trial))
            oracle = community_reference(Matrix(comm), tree)
            for family in (:sorensen, :jaccard)
                full = phylobeta(comm, tree; index_family = family)
                expected = oracle.beta(family)
                for field in (:total, :turnover, :nestedness)
                    @test same_values(getproperty(full, field), getproperty(expected, field))
                end
                pairs = [(rand(rng, 1:128), rand(rng, 1:128)) for _ in 1:100]
                selected = phylobeta(comm, tree; index_family = family, pairs = pairs)
                for field in (:total, :turnover, :nestedness)
                    @test same_values(getproperty(selected, field), [getproperty(full, field)[i, j] for (i, j) in pairs])
                end
            end
            ref = oracle.beta(:jaccard).total
            ref[isnan.(ref)] .= 0
            @test same_values(unifrac(comm, tree), ref)
            ref = 1 .- oracle.beta(:sorensen).total
            ref[isnan.(ref)] .= 0
            @test same_values(phylosor(comm, tree), ref)
        end
    end
end
