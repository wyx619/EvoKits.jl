using Statistics
using CSV

# Independent reference: explicit species lists and BigFloat pair sums,
# without the production traversal, weighted online update or sparse assembly.
function aggregation_trait_reference(matrix, traits, statistic, weighted, corrected)
    n = size(matrix, 1)
    values = fill(NaN, n)
    matched, richness, weights = zeros(Int, n), zeros(Int, n), zeros(n)
    for row in 1:n
        present = findall(>(0), matrix[row, :])
        selected = filter(col -> !ismissing(traits[col]), present)
        richness[row], matched[row] = length(present), length(selected)
        xs = BigFloat[traits[col] for col in selected]
        ws = BigFloat[weighted ? matrix[row, col] : 1 for col in selected]
        total = sum(ws)
        weights[row] = Float64(total)
        if statistic === :sum
            values[row] = Float64(sum(xs .* ws))
        elseif !isempty(selected)
            if statistic === :mean
                values[row] = Float64(sum(xs .* ws) / total)
            elseif !corrected || length(selected) > 1
                pair_values, pair_weights = BigFloat(0), BigFloat(0)
                for i in eachindex(selected), j in (i + 1):length(selected)
                    pair_values += ws[i] * ws[j] * (xs[i] - xs[j])^2
                    pair_weights += ws[i] * ws[j]
                end
                denominator = corrected ? 2pair_weights : total^2
                values[row] = Float64(sqrt(pair_values / denominator))
            end
        end
    end
    return (; values, matched, richness, weights)
end

function aggregation_groups_reference(matrix, assignments, labels, weighted, normalize)
    output = zeros(Float64, size(matrix, 1), length(labels))
    for row in axes(matrix, 1), (index, group) in enumerate(labels)
        selected = [col for col in axes(matrix, 2) if !ismissing(assignments[col]) &&
            string(assignments[col]) == group && matrix[row, col] > 0]
        output[row, index] = sum(weighted ? matrix[row, col] : 1 for col in selected; init = 0.0)
    end
    if normalize
        for row in axes(matrix, 1)
            total = sum(output[row, :])
            total > 0 && (output[row, :] ./= total)
        end
    end
    return output
end

@testset "Community trait summaries: labels, weighting and missing mappings" begin
    matrix = [2 1 0 0; 0 3 1 0; 0 0 0 0; 1 0 0 0]
    species, samples = ["A", "B", "C", "unused"], ["AB", "BC", "empty", "A"]
    traits = Union{Missing, Float64}[2, 6, -2, missing]
    table = DataFrame(taxon = ["C", "extra", "A", "B"], size = [-2, 100, 2, 6])
    mapping = Dict("A" => 2, "B" => 6, "C" => -2, "extra" => 100)
    for input in (matrix, sparse(matrix)), weighted in (false, true), statistic in (:sum, :mean, :sd), corrected in (false, true)
        expected = aggregation_trait_reference(matrix, traits, statistic, weighted, corrected)
        options = (; statistic, abundance_weighted = weighted, corrected)
        data = CommunityData(input; species_labels = species, sample_labels = samples)
        for actual in (map_trait(input, traits; species_labels = species, sample_labels = samples, options...),
            map_trait(data, mapping; options...), map_trait(data, table; species = :taxon, trait = :size, options...))
            @test actual isa CommunityTraitResult
            @test same_values(actual.values, expected.values)
            @test actual.sample_labels == samples
            @test actual.richness == expected.richness
            @test actual.matched_richness == expected.matched
            @test actual.weight_sum == expected.weights
            @test actual.statistic === statistic
            @test actual.abundance_weighted == weighted && actual.corrected == corrected
        end
    end
    partial = Dict("A" => 2, "C" => -2)
    for input in (matrix, sparse(matrix))
        @test_throws ArgumentError map_trait(input, partial; species_labels = species)
        for statistic in (:sum, :mean, :sd)
            actual = map_trait(input, partial; species_labels = species, statistic, missing_traits = :skip)
            expected = aggregation_trait_reference(matrix, [2, missing, -2, missing], statistic, false, true)
            @test same_values(actual.values, expected.values)
            @test actual.richness == [2, 2, 0, 1]
            @test actual.matched_richness == [1, 1, 0, 1]
        end
        @test map_trait(input, Dict(); missing_traits = :skip).values == zeros(4)
        @test all(isnan, map_trait(input, Dict(); statistic = :mean, missing_traits = :skip).values)
    end
    for statistic in (:sum, :mean, :sd)
        empty = map_trait(spzeros(3, 0), Float64[]; statistic)
        @test empty.richness == empty.matched_richness == [0, 0, 0]
        @test empty.weight_sum == [0, 0, 0]
        @test statistic === :sum ? empty.values == [0, 0, 0] : all(isnan, empty.values)
        @test isempty(map_trait(spzeros(0, 2), [2, 4]; statistic).values)
    end
    @test map_trait([1 1 1], [1e16, 1.0, -1e16]).values == [1.0]
    @test map_trait([1 1], [-1e308, 1e308]; statistic = :mean).values == [0.0]
    @test map_trait(reshape([1.0, 1e20], 1, :), [1, 2]; statistic = :sd, abundance_weighted = true).values[1] ≈ sqrt(0.5)
    @test map_trait(ones(1, 1), [1e308]; statistic = :sd, corrected = false).values == [0.0]
    for weighted in (false, true), corrected in (false, true)
        shifted = [1e15, 1e15 + 0.25, 1e15 + 0.5]
        expected = aggregation_trait_reference([1 2 3], shifted, :sd, weighted, corrected)
        @test same_values(map_trait([1 2 3], shifted; statistic = :sd,
            abundance_weighted = weighted, corrected).values, expected.values)
    end
    # Own result arrays and preserve caller data and mappings.
    data = CommunityData(sparse(matrix); species_labels = species, sample_labels = samples)
    saved = deepcopy(data)
    result = map_trait(data, mapping)
    result.values[1] = -99
    result.sample_labels[1] = "changed"
    @test map_trait(data, mapping).values == [8, 4, 0, 2]
    @test data.matrix == saved.matrix && data.sample_labels == saved.sample_labels
    @test mapping == Dict("A" => 2, "B" => 6, "C" => -2, "extra" => 100)
end

@testset "Sparse species groups: counts, abundances, proportions and snapshots" begin
    matrix = [2 1 0 0; 0 3 1 0; 0 0 0 0; 1 0 0 0]
    species, samples = ["A", "B", "C", "unused"], ["AB", "BC", "empty", "A"]
    assignments, labels = ["x", "x", "y", "z"], ["y", "x", "empty", "z"]
    mapping = Dict(zip(species, assignments))
    mapping["extra"] = "outside"
    table = DataFrame(taxon = reverse(species), guild = reverse(assignments))
    for input in (matrix, sparse(matrix)), weighted in (false, true), normalize in (false, true)
        options = (; abundance_weighted = weighted, normalize, group_labels = labels)
        data = CommunityData(input; species_labels = species, sample_labels = samples)
        expected = aggregation_groups_reference(matrix, assignments, labels, weighted, normalize)
        for actual in (aggregate_species(input, assignments; species_labels = species, sample_labels = samples, options...),
            aggregate_species(data, mapping; options...), aggregate_species(data, table; species = :taxon, group = :guild, options...))
            @test actual isa CommunityData
            @test actual.matrix isa SparseMatrixCSC
            @test Matrix(actual.matrix) ≈ expected
            @test actual.sample_labels == samples
            @test actual.species_labels == labels
            @test eltype(actual.matrix) === (normalize ? Float64 : Int)
            @test actual.matrix[3, :] == zeros(4)
        end
    end
    @test aggregate_species(matrix, assignments).species_labels == ["x", "y", "z"]
    @test aggregate_species(matrix, [1, 1, 2, 3]).species_labels == ["1", "2", "3"]
    partial = Dict("A" => "x", "C" => "y")
    @test_throws ArgumentError aggregate_species(matrix, partial; species_labels = species)
    actual = aggregate_species(matrix, partial; species_labels = species, missing_groups = :skip, normalize = true)
    @test Matrix(actual.matrix) == [1 0; 0 1; 0 0; 1 0]
    @test aggregate_species(matrix, Dict(); missing_groups = :skip).matrix == spzeros(Int, 4, 0)
    @test size(aggregate_species(spzeros(0, 2), ["x", "y"])) == (0, 2)
    @test size(aggregate_species(spzeros(3, 0), String[]; group_labels = ["empty"])) == (3, 1)
    @test eltype(aggregate_species(sparse(Bool[1 1]), [1, 1]; abundance_weighted = true).matrix) === Int
    @test aggregate_species(sparse(Float32[1.5 2.5]), [1, 1]; abundance_weighted = true).matrix[1, 1] === 4f0
    @test aggregate_species(sparse([1e308 1e308]), [1, 2]; abundance_weighted = true, normalize = true).matrix == sparse([0.5 0.5])
    explicit = sparse([1, 2, 1], [1, 2, 3], [2.0, 0, 1], 3, 4)
    @test nnz(aggregate_species(explicit, [1, 2, 1, 3]).matrix) == 1
    @test map_trait(explicit, [2, 99, 3, 7]).values == [5, 0, 0]
    data = CommunityData(sparse(matrix); species_labels = species, sample_labels = samples)
    grouped = aggregate_species(data, mapping)
    grouped.matrix[1, 1] = 999
    grouped.species_labels[1] = "changed"
    grouped.sample_labels[1] = "changed"
    @test aggregate_species(data, mapping).matrix[1, 1] == 2
    @test data.species_labels == species && data.sample_labels == samples && data.matrix == matrix
    # Grouped data is usable by existing taxonomic metrics without conversion.
    grouped = aggregate_species(data, mapping)
    @test beta_diss(grouped).total ≈ beta_diss(grouped.matrix).total nans = true
    @test weighted_endemism(grouped).we == weighted_endemism(grouped.matrix).we
end

@testset "Community aggregation validation and finite arithmetic" begin
    for f in (map_trait, aggregate_species)
        @test_throws ArgumentError f(ones(1, 2), [1])
        @test_throws ArgumentError f([-1 1], [1, 2])
        @test_throws ArgumentError f(sparse([NaN 1]), [1, 2])
        @test_throws ArgumentError f(ones(1, 2), [1, 2]; species_labels = ["A", "A"])
        @test_throws ArgumentError f(ones(1, 2), [1, 2]; sample_labels = ["a", "b"])
        @test_throws ArgumentError f(ones(1, 2), (1, 2))
        data = CommunityData(ones(1, 2))
        for keyword in ((species_labels = ["A", "B"],), (sample_labels = ["a"],), (tip_labels = ["A", "B"],))
            @test_throws ArgumentError f(data, [1, 2]; keyword...)
        end
        # Integer and string species identifiers may collide after normalization.
        @test_throws ArgumentError f(ones(1, 1), Dict(1 => 1, "1" => 2))
        @test_throws ArgumentError f(ones(1, 1), DataFrame(species = ["A", "A"], trait = [1, 2], group = [1, 2]))
        @test_throws ArgumentError f(ones(1, 1), DataFrame(wrong = [1]))
    end
    @test_throws ArgumentError map_trait(ones(1, 2), [1, 2]; statistic = :median)
    @test_throws ArgumentError map_trait(ones(1, 2), [1, 2]; missing_traits = :zero)
    for traits in ([1, NaN], [1, Inf], [1, "bad"], [1, nothing], [1, big"1e400"])
        @test_throws ArgumentError map_trait(ones(1, 2), traits)
    end
    @test_throws ArgumentError map_trait([1e308 1e308], [1, 2]; abundance_weighted = true)
    @test_throws ArgumentError map_trait([1 1], [1e308, 1e308])
    @test_throws ArgumentError map_trait(ones(1, 2), [-1e308, 1e308]; statistic = :sd)
    @test_throws ArgumentError map_trait(reshape([big"1e-400"], 1, 1), [1]; abundance_weighted = true)
    @test_throws ArgumentError aggregate_species(ones(1, 2), [1, 2]; missing_groups = :zero)
    @test_throws ArgumentError aggregate_species(ones(1, 2), [1, 2]; group_labels = ["1", "1"])
    @test_throws ArgumentError aggregate_species(ones(1, 2), [1, 2]; group_labels = ["1"])
    for groups in ([1, ""], [1, NaN], [1, 2.0], [1, nothing], Any[1, true])
        @test_throws ArgumentError aggregate_species(ones(1, 2), groups)
    end
    @test_throws ArgumentError aggregate_species([typemax(Int) 1], [1, 1]; abundance_weighted = true)
    @test_throws ArgumentError aggregate_species([1e308 1e308], [1, 1]; abundance_weighted = true)
    @test aggregate_species(ones(1, 2), [1, missing]; missing_groups = :skip).matrix == sparse(reshape([1], 1, 1))
end

@testset "Community aggregation randomized independent formulas and permutations" begin
    rng = MersenneTwister(20261102)
    for trial in 1:35
        n, p = rand(rng, 1:20), rand(rng, 2:40)
        matrix = rand(rng, 0:8, n, p)
        matrix[rand(rng, n, p) .< 0.75] .= 0
        matrix[end, :] .= 0
        species, samples = ["sp_$i" for i in 1:p], ["site_$i" for i in 1:n]
        traits = Union{Missing, Float64}[rand(rng) < 0.2 ? missing : rand(rng, -20:30) / 4 for _ in 1:p]
        groups = Union{Missing, String}[rand(rng) < 0.2 ? missing : "g$(rand(rng, 1:5))" for _ in 1:p]
        labels = ["g$i" for i in 1:6]
        permutation = randperm(rng, p)
        mapping = Dict(species .=> traits)
        grouping = Dict(species .=> groups)
        reordered = CommunityData(sparse(matrix[:, permutation]); species_labels = species[permutation], sample_labels = samples)
        saved = copy(reordered.matrix)
        for weighted in (false, true)
            for statistic in (:sum, :mean, :sd), corrected in (false, true)
                expected = aggregation_trait_reference(matrix, traits, statistic, weighted, corrected)
                options = (; statistic, abundance_weighted = weighted, corrected, missing_traits = :skip)
                for actual in (map_trait(matrix, traits; options...), map_trait(reordered, mapping; options...))
                    @test same_values(actual.values, expected.values)
                    @test actual.richness == expected.richness
                    @test actual.matched_richness == expected.matched
                    @test actual.weight_sum == expected.weights
                end
                if weighted && statistic !== :sum
                    scaled = map_trait(CommunityData(reordered.matrix .* 3; species_labels = reordered.species_labels), mapping; options...)
                    @test same_values(scaled.values, expected.values)
                end
            end
            for normalize in (false, true)
                expected = aggregation_groups_reference(matrix, groups, labels, weighted, normalize)
                for actual in (aggregate_species(matrix, groups; group_labels = labels, missing_groups = :skip, abundance_weighted = weighted, normalize),
                    aggregate_species(reordered, grouping; group_labels = labels, missing_groups = :skip, abundance_weighted = weighted, normalize))
                    @test Matrix(actual.matrix) ≈ expected
                    @test actual.species_labels == labels
                    @test actual.matrix isa SparseMatrixCSC
                    @test eltype(actual.matrix) === (normalize ? Float64 : Int)
                end
            end
        end
        @test reordered.matrix == saved
    end
end

@testset "Community aggregation saved local R reference" begin
    inputs = CSV.read(joinpath(@__DIR__, "fixtures", "community_aggregation_input.tsv"), DataFrame; delim = '\t')
    trait_reference = CSV.read(joinpath(@__DIR__, "fixtures", "community_aggregation_traits_r.tsv"), DataFrame;
        delim = '\t', missingstring = "NaN")
    group_reference = CSV.read(joinpath(@__DIR__, "fixtures", "community_aggregation_groups_r.tsv"), DataFrame;
        delim = '\t', missingstring = "NaN")
    for case in unique(inputs.case)
        table = inputs[inputs.case .== case, :]
        traits = unique(table[:, [:species, :trait]])
        groups = unique(table[:, [:species, :group]])
        rtraits = trait_reference[trait_reference.case .== case, :]
        rgroups = group_reference[group_reference.case .== case, :]
        data = CommunityData(table; sample = :site, sample_labels = rtraits.site)
        labels = unique(rgroups.group)
        for matrix in (Matrix(data.matrix), data.matrix)
            community = CommunityData(matrix; sample_labels = data.sample_labels, species_labels = data.species_labels)
            for (column, statistic, weighted) in ((:presence_sum, :sum, false), (:presence_mean, :mean, false),
                (:presence_sd, :sd, false), (:weighted_sum, :sum, true), (:weighted_mean, :mean, true))
                result = map_trait(community, traits; statistic, abundance_weighted = weighted)
                reference = Float64.(coalesce.(rtraits[!, column], NaN))
                # R drops empty sites. Their explicit retention is verified
                # independently above; R comparisons cover nonempty sites.
                mask = rtraits.richness .> 0
                @test same_values(result.values[mask], reference[mask])
                @test result.richness == rtraits.richness
            end
            for (column, weighted, normalize) in ((:count, false, false), (:abundance, true, false),
                (:count_proportion, false, true), (:abundance_proportion, true, true))
                result = aggregate_species(community, groups; group_labels = labels,
                    abundance_weighted = weighted, normalize)
                reference = reshape(Float64.(coalesce.(rgroups[!, column], NaN)), size(result))
                reference[isnan.(reference)] .= 0 # new empty-row normalization convention
                @test same_values(Matrix(result.matrix), reference)
                @test result.sample_labels == rtraits.site && result.species_labels == labels
            end
        end
    end
end

@testset "Community aggregation large sparse input does not densify" begin
    rng = MersenneTwister(420)
    n, p = 2000, 8000
    matrix = sparse(rand(rng, 1:n, 6000), rand(rng, 1:p, 6000), ones(Int, 6000), n, p)
    traits, groups = rand(rng, p), [mod1(i, 8) for i in 1:p]
    for statistic in (:sum, :mean, :sd)
        map_trait(matrix, traits; statistic, abundance_weighted = true)
        @test (@allocated map_trait(matrix, traits; statistic, abundance_weighted = true)) < 10_000_000
    end
    aggregate_species(matrix, groups; abundance_weighted = true, normalize = true)
    @test (@allocated aggregate_species(matrix, groups; abundance_weighted = true, normalize = true)) < 10_000_000
    grouped = aggregate_species(matrix, groups; abundance_weighted = true)
    @test size(grouped) == (n, 8)
    @test nnz(grouped.matrix) <= nnz(matrix)
    @test sum(grouped.matrix) == sum(matrix)
    @test Base.summarysize(grouped) < 1_000_000
end
