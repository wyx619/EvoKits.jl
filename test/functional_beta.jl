using Test, EvoKits, DataFrames, Random, SparseArrays, CSV, StableRNGs

# Recompute the objective from original traits and independently fitted means
# and modes, without using the implementation's stored prototypes.
function functional_test_objective(numeric, categorical, membership; numeric_weights, gamma, categorical_weights)
    result = BigFloat(0)
    for group in unique(membership)
        members = findall(==(group), membership)
        for j in axes(numeric, 2)
            values = BigFloat.(numeric[members, j])
            center = sum(values) / length(values)
            result += BigFloat(numeric_weights[j]) * sum((x - center)^2 for x in values)
        end
        for j in axes(categorical, 2)
            values = categorical[members, j]
            largest = maximum(count(==(value), values) for value in unique(values))
            result += BigFloat(gamma) * BigFloat(categorical_weights[j]) * (length(values) - largest)
        end
    end
    return Float64(result)
end

@testset "Functional grouping: mixed traits and independent objective" begin
    traits = DataFrame(species = ["A", "B", "C", "D"],
        height = [0.0, 2.0, 10.0, 12.0], form = ["a", "a", "b", "b"])
    original = deepcopy(traits)
    grouping = functional_groups(traits; k = 2, numeric_scale = :none,
        gamma = 10, nstart = 20, rng = MersenneTwister(42))
    @test grouping isa FunctionalGroupingResult
    @test grouping.membership == [1, 1, 2, 2]
    @test grouping.group_sizes == [2, 2]
    @test grouping.group_labels == ["functional_1", "functional_2"]
    @test grouping.numeric_traits == ["height"] && grouping.categorical_traits == ["form"]
    @test grouping.numeric_centers[:, 1] == [1, 11]
    @test grouping.categorical_modes[:, 1] == ["a", "b"]
    @test grouping.objective == 4
    @test grouping.converged && grouping.iterations >= 1 && grouping.nstart == 20
    @test isequal(traits, original)
    repeat = functional_groups(traits; k = 2, numeric_scale = :none,
        gamma = 10, nstart = 20, rng = MersenneTwister(42))
    @test repeat.membership == grouping.membership && repeat.objective == grouping.objective

    # Exhaustively enumerate all nonempty two-group partitions of this small
    # dataset; this verifies the reported minimum for a known separated case.
    minimum_objective = Inf
    for bits in 1:(2^4 - 2)
        partition = [iszero(bits & (1 << (i - 1))) ? 1 : 2 for i in 1:4]
        objective = functional_test_objective(reshape(traits.height, :, 1),
            reshape(traits.form, :, 1), partition;
            numeric_weights = [1], gamma = 10, categorical_weights = [1])
        minimum_objective = min(minimum_objective, objective)
    end
    @test grouping.objective == minimum_objective

    numeric = DataFrame(species = ["x", "y", "z"], x = [2.0, 3.0, 9.0])
    only_numeric = functional_groups(numeric; k = 1, numeric_scale = :none,
        gamma = :auto, rng = MersenneTwister(1))
    @test only_numeric.gamma == 1 && only_numeric.categorical_traits == String[]
    @test only_numeric.objective ≈ 86 / 3
    categorical = DataFrame(species = ["x", "y", "z", "w"], type = ["b", "a", "b", "a"])
    only_categorical = functional_groups(categorical; k = 1, gamma = :auto,
        rng = MersenneTwister(1))
    @test only_categorical.gamma == 1 && only_categorical.numeric_traits == String[]
    @test only_categorical.objective == 2 && only_categorical.categorical_modes[1, 1] == "a"
    tied = DataFrame(species = ["x", "y", "z", "w"],
        x = fill(5.0, 4), type = fill("same", 4))
    groups = functional_groups(tied; k = 4, rng = MersenneTwister(2))
    @test groups.converged && groups.group_sizes == ones(Int, 4)
    @test groups.objective == 0 && groups.gamma == 1
    @test groups.numeric_scales == [1]
    @test functional_groups(tied; k = 2, nstart = 1, maxiter = 1,
        rng = MersenneTwister(2)).converged == false

    coded = DataFrame(species = ["a", "b", "c"], code = [1, 1, 2], alive = [true, false, true])
    codes = functional_groups(coded; k = 1, categorical_traits = (:code,), gamma = 2,
        trait_weights = Dict(:code => 3, "alive" => 0.5), rng = MersenneTwister(2))
    @test codes.numeric_traits == String[] && codes.categorical_traits == ["code", "alive"]
    @test codes.categorical_weights == [3, 0.5]
    @test codes.objective == 7
    @test codes.categorical_modes == reshape(["1", "true"], 1, 2)

    # Mean pairwise numeric and categorical contributions agree under :auto.
    auto = functional_groups(traits; k = 2, numeric_scale = :none,
        gamma = :auto, nstart = 2, rng = MersenneTwister(4))
    numeric_pairs = sum((traits.height[i] - traits.height[j])^2 for i in 1:3 for j in (i + 1):4)
    categorical_pairs = sum(traits.form[i] != traits.form[j] for i in 1:3 for j in (i + 1):4)
    @test auto.gamma ≈ numeric_pairs / categorical_pairs
end

@testset "Functional grouping: randomized objective and reproducibility" begin
    rng = MersenneTwister(71)
    for trial in 1:35
        n = rand(rng, 4:10)
        numeric = rand(rng, -8:8, n, 2)
        categories = [string(rand(rng, 'a':'c')) for _ in 1:n, _ in 1:2]
        traits = DataFrame(species = ["sp$i" for i in 1:n],
            x = numeric[:, 1], y = numeric[:, 2], a = categories[:, 1], b = categories[:, 2])
        k = rand(rng, 1:min(4, n))
        options = (; k, numeric_scale = :none, gamma = 1.5,
            trait_weights = Dict(:x => 2.0, :y => 0.5, :a => 3.0, :b => 2.0), maxiter = 30)
        one = functional_groups(traits; options..., nstart = 1, rng = MersenneTwister(trial))
        result = functional_groups(traits; options..., nstart = 5, rng = MersenneTwister(trial))
        expected = functional_test_objective(numeric, categories, result.membership;
            numeric_weights = [2, 0.5], gamma = 1.5, categorical_weights = [3, 2])
        @test result.objective ≈ expected atol = 1e-10
        @test result.objective <= one.objective + 1e-10
        @test sum(result.group_sizes) == n && all(>(0), result.group_sizes)
        @test result.group_sizes == [count(==(i), result.membership) for i in 1:k]
        @test result.numeric_offsets == Float64.(numeric[1, :])
        for group in 1:k
            members = findall(==(group), result.membership)
            for j in 1:2
                expected_center = sum(numeric[members, j] .- result.numeric_offsets[j]) / length(members)
                @test result.numeric_centers[group, j] ≈ expected_center atol = 1e-12
                modes = categories[members, j]
                frequencies = Dict(value => count(==(value), modes) for value in unique(modes))
                largest = maximum(values(frequencies))
                @test result.categorical_modes[group, j] == minimum([x for (x, v) in frequencies if v == largest])
            end
        end
    end
    large_offset = 1e20
    trait = DataFrame(species = ["a", "b", "c"],
        val = [large_offset, nextfloat(large_offset), nextfloat(nextfloat(large_offset))])
    fitted = functional_groups(trait; k = 1, rng = MersenneTwister(1))
    @test fitted.numeric_scales[1] > 0 && isfinite(fitted.objective)
    @test fitted.objective ≈ 2
end

@testset "Functional beta: sparse aggregation, reuse and integration" begin
    traits = DataFrame(species = ["C", "A", "D", "B", "extra"],
        height = [10.0, 0.0, 12.0, 2.0, 1000.0], form = ["b", "a", "b", "a", "z"])
    matrix = sparse([2 3 0 0; 0 0 1 4; 1 0 1 0; 0 0 0 0])
    data = CommunityData(matrix; species_labels = ["A", "B", "C", "D"],
        sample_labels = ["west", "east", "both", "empty"])
    group = functional_groups(data, traits; k = 2, numeric_scale = :none,
        gamma = 10, nstart = 12, rng = MersenneTwister(42))
    @test group.species_labels == data.species_labels
    @test group.membership == [1, 1, 2, 2]
    @test Matrix(collapse_traits(data, group).matrix) == [1 0; 0 1; 1 1; 0 0]
    @test Matrix(collapse_traits(data, group; statistic = :richness).matrix) == [2 0; 0 2; 1 1; 0 0]
    @test Matrix(collapse_traits(data, group; statistic = :abundance).matrix) == [5 0; 0 5; 1 1; 0 0]
    @test Matrix(collapse_traits(data, group; normalize = true).matrix) == [1 0; 0 1; 0.5 0.5; 0 0]
    @test collapse_traits(matrix, group).matrix == collapse_traits(data, group).matrix
    result = functional_beta(data, group)
    @test result isa FunctionalBetaResult
    @test result.community.matrix isa SparseMatrixCSC
    @test result.beta.sample_labels == data.sample_labels
    @test isempty(result.excluded_species)
    @test result.beta.total[1, 2] == 1 && result.beta.turnover[1, 2] == 1
    @test result.beta.total[1, 3] ≈ 1 / 3 && result.beta.turnover[1, 3] == 0
    @test isnan(result.beta.total[4, 4])
    direct = functional_beta(data, traits; k = 2, numeric_scale = :none,
        gamma = 10, nstart = 12, rng = MersenneTwister(42))
    @test isequal(direct.beta.total, result.beta.total)
    @test direct.grouping.membership == group.membership
    abundance = functional_beta(data, group; statistic = :abundance, index_family = :jaccard,
        pairs = [(1, 2), (1, 3), (4, 4)])
    @test Matrix(abundance.community.matrix) == [5 0; 0 5; 1 1; 0 0]
    @test isequal(abundance.beta.total, [1.0, 0.5, NaN])
    @test abundance.beta.pairs == [(1, 2), (1, 3), (4, 4)]
    matrix_result = functional_beta(matrix, group; species_labels = data.species_labels,
        sample_labels = data.sample_labels)
    @test isequal(matrix_result.beta.total, result.beta.total)
    @test isequal(functional_beta(matrix, group).beta.total, result.beta.total)
    subset = CommunityData(matrix[[1, 2, 3], :]; species_labels = data.species_labels,
        sample_labels = data.sample_labels[[1, 2, 3]])
    @test length(phyloregion(functional_beta(subset, group).beta; k = 2).membership) == 3
    shuffled = CommunityData(matrix[:, [4, 2, 1, 3]];
        species_labels = data.species_labels[[4, 2, 1, 3]], sample_labels = data.sample_labels)
    @test collapse_traits(shuffled, group).matrix == result.community.matrix
    result.grouping.membership[1] = 2
    @test group.membership == [1, 1, 2, 2]
    result.community.matrix[1, 1] = 0
    @test data.matrix == matrix && collapse_traits(data, group).matrix[1, 1] == 1
end

@testset "Functional grouping: missing values, validation and ownership" begin
    traits = DataFrame(species = ["a", "b", "c", "extra"],
        value = Union{Missing, Float64}[1, missing, 3, 4], form = ["x", "x", "y", "z"])
    data = CommunityData(sparse([1 0 1 0; 0 0 0 0]);
        species_labels = ["a", "b", "c", "missing"], sample_labels = ["s1", "empty"])
    @test_throws ArgumentError functional_groups(data, traits; k = 2)
    group = functional_groups(data, traits; k = 2, missing_traits = :skip, rng = MersenneTwister(3))
    @test group.species_labels == ["a", "c"]
    @test group.excluded_species == ["b", "missing"]
    @test_throws ArgumentError collapse_traits(data, group)
    @test Matrix(collapse_traits(data, group; missing_traits = :skip).matrix) == [1 1; 0 0]
    fit = functional_beta(data, traits; k = 2, missing_traits = :skip, rng = MersenneTwister(3))
    @test fit.grouping.excluded_species == ["b", "missing"]
    @test fit.excluded_species == ["b", "missing"]
    @test functional_beta(data, group; missing_traits = :skip).excluded_species == ["b", "missing"]
    @test_throws ArgumentError functional_groups(traits; k = 2)
    @test_throws ArgumentError functional_groups(data, traits; k = 3, missing_traits = :skip)
    complete = dropmissing(traits)
    for k in (0, -1, true, 4, 2.5)
        @test_throws ArgumentError functional_groups(complete; k)
    end
    for gamma in (0, -1, Inf, NaN, true, :bad)
        @test_throws ArgumentError functional_groups(complete; k = 1, gamma)
    end
    @test_throws ArgumentError functional_groups(complete; k = 1, nstart = 0)
    @test_throws ArgumentError functional_groups(complete; k = 1, maxiter = 0)
    @test_throws ArgumentError functional_groups(complete; k = 1, numeric_scale = :bad)
    @test_throws ArgumentError functional_groups(complete; k = 1, missing_traits = :bad)
    @test_throws ArgumentError functional_groups(complete; k = 1, trait_weights = Dict(:nope => 1))
    @test_throws ArgumentError functional_groups(complete; k = 1, trait_weights = Dict(:value => 0))
    @test_throws ArgumentError functional_groups(complete; k = 1, trait_cols = [:value, :value])
    @test_throws ArgumentError functional_groups(complete; k = 1, categorical_traits = [:nope])
    @test_throws ArgumentError functional_groups(DataFrame(species = ["a", "a"], x = [1, 2]); k = 1)
    @test_throws ArgumentError functional_groups(DataFrame(species = ["a"], x = [Inf]); k = 1)
    @test_throws ArgumentError functional_groups(DataFrame(species = ["a"], x = [""]); k = 1)
    @test_throws ArgumentError functional_groups(DataFrame(species = ["a", "b"], x = Any[1, "1"]); k = 1)
    @test_throws ArgumentError functional_groups(DataFrame(species = ["a", "b"], x = Any[1, 1.0]);
        k = 1, categorical_traits = (:x,))
    @test_throws ArgumentError collapse_traits(data, group; statistic = :bad, missing_traits = :skip)
    damaged = deepcopy(group)
    damaged.membership[1] = 0
    @test_throws ArgumentError collapse_traits(data, damaged; missing_traits = :skip)
    damaged = deepcopy(group)
    damaged.group_sizes[1] += 1
    @test_throws ArgumentError collapse_traits(data, damaged; missing_traits = :skip)
end

@testset "Functional grouping: automatic gamma avoids species-pair storage" begin
    n = 3000
    traits = DataFrame(species = ["sp$i" for i in 1:n], x = collect(1:n),
        type = [isodd(i) ? "a" : "b" for i in 1:n])
    functional_groups(traits; k = 3, nstart = 1, maxiter = 2, rng = MersenneTwister(1))
    bytes = @allocated functional_groups(traits; k = 3, nstart = 1, maxiter = 2, rng = MersenneTwister(1))
    @test bytes < 15_000_000
end

@testset "Functional beta: saved local clustMixType and phyloregion references" begin
    path = joinpath(@__DIR__, "fixtures", "functional")
    saved_groups = CSV.read(joinpath(path, "kproto_r.tsv"), DataFrame; delim = '\t')
    saved_community = CSV.read(joinpath(path, "community_r.tsv"), DataFrame; delim = '\t')
    saved_beta = CSV.read(joinpath(path, "beta_r.tsv"), DataFrame; delim = '\t')
    for id in unique(saved_groups.case)
        rows = saved_groups[saved_groups.case .== id, :]
        cfg = first(eachrow(rows))
        traits = select(rows, :species, :x, :y, :a, :b)
        comm = saved_community[saved_community.case .== id, :]
        matrix = permutedims(reshape(comm.value, nrow(rows), 7))
        data = CommunityData(sparse(matrix); species_labels = rows.species,
            sample_labels = ["site_$i" for i in 1:7])
        result = functional_groups(data, traits; k = cfg.k, numeric_scale = Symbol(cfg.scale),
            gamma = cfg.gamma_mode == "auto" ? :auto : cfg.gamma,
            trait_weights = Dict(:x => cfg.wx, :y => cfg.wy, :a => cfg.wa, :b => cfg.wb),
            nstart = 1, maxiter = 200, rng = StableRNG(cfg.seed))
        @test result.membership == rows.group
        @test result.species_labels == rows.species
        @test result.group_sizes == [count(==(g), rows.group) for g in 1:cfg.k]
        @test result.objective ≈ cfg.objective atol = 1e-10 rtol = 1e-11
        @test result.gamma ≈ cfg.gamma rtol = 1e-12
        @test result.numeric_offsets ≈ [cfg.offset_x, cfg.offset_y] atol = 1e-12
        @test result.numeric_scales ≈ [cfg.scale_x, cfg.scale_y] atol = 1e-12
        @test result.numeric_centers[result.membership, 1] ≈ rows.center_x atol = 1e-11
        @test result.numeric_centers[result.membership, 2] ≈ rows.center_y atol = 1e-11
        @test result.categorical_modes[result.membership, 1] == rows.mode_a
        @test result.categorical_modes[result.membership, 2] == rows.mode_b
        @test result.converged
        # Values repeat per input species so all input zero columns and group
        # memberships remain visible in a single, small community fixture.
        expected = unique(select(comm, :site, :group, :presence, :richness, :abundance, :proportion))
        for statistic in (:presence, :richness, :abundance)
            collapsed = collapse_traits(data, result; statistic)
            @test collapsed.matrix isa SparseMatrixCSC
            @test [collapsed.matrix[row.site, row.group] for row in eachrow(expected)] == expected[!, statistic]
        end
        normalized = collapse_traits(data, result; statistic = :abundance, normalize = true)
        @test [normalized.matrix[row.site, row.group] for row in eachrow(expected)] ≈ expected.proportion atol = 1e-13
        for family in (:sorensen, :jaccard)
            reference = saved_beta[(saved_beta.case .== id) .& (saved_beta.family .== String(family)), :]
            beta = functional_beta(data, result; index_family = family).beta
            for field in (:total, :turnover, :nestedness)
                @test isapprox([getfield(beta, field)[row.i, row.j] for row in eachrow(reference)],
                    reference[!, field]; nans = true, atol = 1e-13)
            end
            pairs = [(row.i, row.j) for row in eachrow(reference)]
            selected = functional_beta(data, result; index_family = family, pairs).beta
            @test isapprox(selected.total, reference.total; nans = true, atol = 1e-13)
        end
    end
end
