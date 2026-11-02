using CSV

# Independent Boolean oracle; nestedness uses its defining fraction rather
# than subtracting the production total and turnover results.
function taxonomic_beta_reference(community, family)
    presence = community .> 0
    n = size(community, 1)
    total, turnover, nestedness = (fill(NaN, n, n) for _ in 1:3)
    for i in 1:n, j in 1:n
        left, right = presence[i, :], presence[j, :]
        a = count(left .& right)
        b = count(left .& .!right)
        c = count(.!left .& right)
        m = min(b, c)
        d = family === :sorensen ? 2a + b + c : a + b + c
        v = family === :sorensen ? a + m : a + 2m
        total[i, j] = d > 0 ? (b + c) / d : NaN
        turnover[i, j] = v > 0 ? (family === :sorensen ? m : 2m) / v : NaN
        nestedness[i, j] = d > 0 && v > 0 ? (abs(b - c) / d) * (a / v) : NaN
    end
    return (; total, turnover, nestedness)
end

@testset "Taxonomic beta: hand calculations and empty sites" begin
    comm = [1 1 1 0; 1 0 0 1; 1 0 0 0; 0 0 0 1; 0 0 0 0]
    labels = ["ABC", "AD", "A", "D", "empty"]
    for input in (comm, Float64.(comm), comm .> 0, sparse(comm), view(comm, :, :))
        sor = beta_diss(input; sample_labels = labels)
        jac = beta_diss(input; index_family = :jaccard)
        @test sor isa BetaDissResult{Matrix{Float64}}
        @test sor.sample_labels == labels
        @test sor.sample_labels !== labels
        @test sor.index_family === :sorensen
        @test sor.pairs === nothing
        @test jac.sample_labels == ["sample_$i" for i in 1:5]
        @test sor.total[1, 2] ≈ 3 / 5
        @test sor.turnover[1, 2] ≈ 1 / 2
        @test sor.nestedness[1, 2] ≈ 1 / 10
        @test jac.total[1, 2] ≈ 3 / 4
        @test jac.turnover[1, 2] ≈ 2 / 3
        @test jac.nestedness[1, 2] ≈ 1 / 12
        # Strict nesting and disjoint nonempty sites.
        @test sor.total[1, 3] == sor.nestedness[1, 3] == 1 / 2
        @test sor.turnover[1, 3] == 0
        @test jac.total[1, 3] ≈ 2 / 3
        @test jac.nestedness[1, 3] ≈ 2 / 3
        for result in (sor, jac)
            @test result.total[1, 4] == result.turnover[1, 4] == 1
            @test result.nestedness[1, 4] == 0
            @test all(iszero, [getproperty(result, f)[i, i] for f in
                (:total, :turnover, :nestedness), i in 1:4])
            @test result.total[1, 5] == 1
            @test isnan(result.turnover[1, 5]) && isnan(result.nestedness[1, 5])
            @test all(isnan, [getproperty(result, f)[5, 5] for f in
                (:total, :turnover, :nestedness)])
        end
        pairs = [(2, 1), (1, 2), (1, 1), (5, 5), (5, 2), (2, 1)]
        for family in (:sorensen, :jaccard)
            full = family === :sorensen ? sor : jac
            selected = beta_diss(input; index_family = family,
                pairs = (p for p in pairs), sample_labels = labels)
            @test selected isa BetaDissResult{Vector{Float64}}
            @test selected.sample_labels == labels
            @test selected.pairs == pairs
            @test selected.pairs !== pairs
            for field in (:total, :turnover, :nestedness)
                @test same_values(getproperty(selected, field),
                    [getproperty(full, field)[i, j] for (i, j) in pairs])
            end
            empty = beta_diss(input; index_family = family, pairs = [])
            @test isempty(empty.total) && isempty(empty.turnover) && isempty(empty.nestedness)
            @test empty.pairs == []
        end
    end
    # Magnitudes never weight presence, including subnormal positives and
    # explicitly stored zeros. Inputs and returned result arrays are independent.
    input = sparse([1, 1, 2, 2, 3], [1, 2, 1, 3, 2],
        [0.0, nextfloat(0.0), 1e308, 2.0, 1e-20], 4, 4)
    before = copy(input)
    actual = beta_diss(input)
    expected = beta_diss(Matrix(input) .> 0)
    for field in (:total, :turnover, :nestedness)
        @test same_values(getproperty(actual, field), getproperty(expected, field))
    end
    @test input == before && nnz(input) == nnz(before)
    @test nonzeros(input) == nonzeros(before)
    actual.total[1, 2] = -1
    @test actual.turnover[1, 2] >= 0
    @test same_values(beta_diss(input).total, expected.total)
    for input in (zeros(0, 3), spzeros(0, 3))
        @test size(beta_diss(input).total) == (0, 0)
    end
    for input in (zeros(3, 0), spzeros(3, 0), spzeros(128, 0), spzeros(128, 4))
        result = beta_diss(input)
        @test all(isnan, result.total)
        @test all(isnan, result.turnover)
        @test all(isnan, result.nestedness)
    end
    @test beta_diss(comm; pairs = [[1, 2], [2, 1]]).pairs == [(1, 2), (2, 1)]
    @test_throws ArgumentError beta_diss(comm; index_family = :wrong)
    @test_throws ArgumentError beta_diss(comm; sample_labels = ["wrong"])
    for bad in ([(0, 1)], [(1, 6)], [(1.0, 2)], [(true, 2)], [(1,)], [1], [(1, 2, 3)])
        @test_throws ArgumentError beta_diss(comm; pairs = bad)
    end
    for bad in (-1.0, Inf, -Inf, NaN), input in ([bad 0], sparse([bad 0]))
        @test_throws ArgumentError beta_diss(input)
        @test_throws ArgumentError beta_diss(input; pairs = [])
    end
end

@testset "Taxonomic beta: independent randomized oracle and permutations" begin
    rng = MersenneTwister(20261103)
    for trial in 1:50
        n = trial <= 4 ? (1, 127, 128, 140)[trial] : rand(rng, 2:20)
        p = rand(rng, 3:40)
        comm = (rand(rng, n, p) .< rand(rng)) .* rand(rng, 1:9, n, p)
        comm[end, :] .= 0
        comm[:, end] .= 0
        if n > 3
            comm[1, 1] = 1
            comm[2, :] .= comm[1, :]
        end
        rowperm, colperm = randperm(rng, n), randperm(rng, p)
        pairs = [(rand(rng, 1:n), rand(rng, 1:n)) for _ in 1:25]
        for family in (:sorensen, :jaccard)
            expected = taxonomic_beta_reference(comm, family)
            for input in (comm, sparse(comm))
                full = beta_diss(input; index_family = family)
                selected = beta_diss(input; index_family = family, pairs = pairs)
                reordered = beta_diss(input[rowperm, colperm]; index_family = family)
                for field in (:total, :turnover, :nestedness)
                    values = getproperty(full, field)
                    oracle = getproperty(expected, field)
                    @test same_values(values, oracle)
                    @test same_values(values, transpose(values))
                    @test all(x -> isnan(x) || 0 <= x <= 1, values)
                    @test same_values(getproperty(selected, field), [oracle[i, j] for (i, j) in pairs])
                    @test same_values(getproperty(reordered, field), oracle[rowperm, rowperm])
                end
                valid = isfinite.(full.nestedness)
                @test full.total[valid] ≈ (full.turnover .+ full.nestedness)[valid]
                @test same_values(full.total, beta_diss(input .> 0; index_family = family).total)
            end
        end
    end
end

@testset "Taxonomic beta equals phylogenetic beta on a unit star tree" begin
    mktempdir() do dir
        path = joinpath(dir, "star.nwk")
        write_tree(path, "(A:1,B:1,C:1,D:1);")
        tree = read_tree(path)
        comm = repeat([1 1 0 0; 0 1 1 1; 0 0 0 0; 1 0 0 0], 32, 1)
        for family in (:sorensen, :jaccard), pairs in (nothing, [(2, 1), (3, 3), (128, 1)])
            tax = beta_diss(sparse(comm); index_family = family, pairs = pairs)
            phy = phylobeta(sparse(comm), tree; index_family = family, pairs = pairs)
            for field in (:total, :turnover, :nestedness)
                @test same_values(getproperty(tax, field), getproperty(phy, field))
            end
        end
    end
end

@testset "Taxonomic selected pairs retain sparse storage" begin
    # A dense Float64 community would occupy 144 MB; three full result matrices
    # would occupy 54 MB. Warm the selected-pair method before measuring it.
    rng = MersenneTwister(20261104)
    comm = sparse(rand(rng, 1:1500, 4000), rand(rng, 1:12000, 4000), ones(4000), 1500, 12000)
    pairs = [(1, 1500), (20, 40), (30, 30)]
    beta_diss(comm; pairs = pairs)
    bytes = @allocated beta_diss(comm; pairs = pairs)
    @test bytes < 8_000_000
end

@testset "Taxonomic beta: saved phyloregion and betapart outputs" begin
    # phyloregion 1.0.9 beta_diss and betapart 1.6.1 beta.pair agree on this
    # fixture. The saved output is read directly; R is not needed to run tests.
    comm = [1 1 1 0 0; 1 0 0 1 0; 1 0 0 0 0;
            0 0 0 1 0; 1 1 1 1 0; 0 1 0 0 1]
    expected = CSV.read(joinpath(@__DIR__, "fixtures", "taxonomic_beta_r.tsv"),
        DataFrame; delim = '\t')
    for family in (:sorensen, :jaccard)
        reference = expected[expected.family .== String(family), :]
        pairs = collect(zip(reference.i, reference.j))
        for input in (comm, sparse(comm), sparse(repeat(comm, 22, 1)))
            full = beta_diss(input; index_family = family)
            selected = beta_diss(input; index_family = family, pairs = pairs)
            for field in (:total, :turnover, :nestedness)
                @test same_values([getproperty(full, field)[i, j] for (i, j) in pairs], reference[!, field])
                @test same_values(getproperty(selected, field), reference[!, field])
            end
        end
    end
end
