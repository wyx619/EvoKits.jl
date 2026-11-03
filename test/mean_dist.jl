using Test, EvoKits, Random, StableRNGs, SparseArrays, DataFrames, CSV

# Reference each requested pair separately in arbitrary precision, without
# using the production alignment, accumulator, symmetry or pair-index helpers.
function mean_dist_reference(inputs, labels, selected; propagate = false)
    means, counts = Float64[], Int[]
    for (a, b) in selected
        observations = BigFloat[]
        undefined = false
        for input in inputs
            i, j = findfirst(==(labels[a]), input.sample_labels),
                findfirst(==(labels[b]), input.sample_labels)
            (i === nothing || j === nothing) && continue
            if haskey(input, :pairs)
                position = findfirst(p -> p == (i, j) || p == (j, i), input.pairs)
                position === nothing && continue
                value = input.distances[position]
            else
                value = input.distances[i, j]
            end
            if isnan(value)
                undefined = true
            else
                push!(observations, BigFloat(value))
            end
        end
        push!(counts, length(observations))
        push!(means, isempty(observations) || (propagate && undefined) ? NaN :
            Float64(sum(observations) / length(observations)))
    end
    return means, counts
end

function mean_dist_test_tree(newick)
    return mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, newick)
        read_tree(path)
    end
end

@testset "Mean distances: label alignment, subsets and undefined observations" begin
    labels = ["c", "a", "b", "never"]
    first_input = (distances = [0.0 2 4; 2 0 NaN; 4 NaN NaN],
        sample_labels = ["c", "a", "b"])
    second = (distances = [0.0 8; 8 NaN], sample_labels = ["b", "a"])
    third = (distances = [0.0 NaN; NaN 0], sample_labels = ["a", "c"])
    inputs = [first_input, second, third]
    average = mean_dist(inputs; sample_labels = labels)
    @test average isa MeanDistResult{2}
    @test average.sample_labels == labels
    @test average.component === :distance
    @test average.index_family === nothing
    @test average.replicates == 3
    @test average.pairs === nothing
    @test isequal(average.distances, [0.0 2 4 NaN; 2 0 8 NaN; 4 8 0 NaN; NaN NaN NaN NaN])
    @test average.counts == [2 1 1 0; 1 2 1 0; 1 1 1 0; 0 0 0 0]
    propagated = mean_dist(inputs; sample_labels = labels, nan_policy = :propagate)
    @test propagated.counts == average.counts
    @test isnan(propagated.distances[1, 2])
    @test isnan(propagated.distances[2, 3])
    @test isnan(propagated.distances[2, 2])
    @test isnan(propagated.distances[3, 3])
    @test propagated.distances[1, 3] == 4
    @test mean_dist(inputs).sample_labels == first_input.sample_labels
    @test isequal(mean_dist(reverse(inputs); sample_labels = labels).distances, average.distances)
    @test mean_dist([(distances = zeros(0, 0), sample_labels = String[])]).counts == zeros(Int, 0, 0)
    empty_input = mean_dist([(distances = zeros(0, 0), sample_labels = String[])]; sample_labels = labels)
    @test all(isnan, empty_input.distances)
    @test all(iszero, empty_input.counts)
    @test empty_input.replicates == 1
    integer_labels = mean_dist([(distances = [0 3; 3 0], sample_labels = [7, 2])]; sample_labels = [2, 7])
    @test integer_labels.sample_labels == ["2", "7"]
    @test integer_labels.distances == [0 3; 3 0]
    @test mean_dist([sparse([0 2; 2 0])]).distances == [0 2; 2 0]
end

@testset "Mean distances: pair vectors, repetitions and partial coverage" begin
    labels = ["c", "a", "b", "never"]
    full = (distances = [0.0 2 4; 2 0 NaN; 4 NaN NaN], sample_labels = labels[1:3])
    # Source pair indices refer to its own reordered labels. (b,a) occurs twice.
    sparse_pairs = (distances = [8.0, 8, 0], sample_labels = ["b", "a"],
        pairs = [(1, 2), (2, 1), (1, 1)])
    requested = [(2, 3), (3, 2), (1, 1), (4, 1), (2, 2), (3, 3), (2, 3)]
    average = mean_dist([full, sparse_pairs]; sample_labels = labels, pairs = requested)
    @test average isa MeanDistResult{1}
    @test average.pairs == requested
    @test isequal(average.distances, [8.0, 8, 0, NaN, 0, 0, 8])
    @test average.counts == [1, 1, 1, 0, 1, 1, 1]
    expected, counts = mean_dist_reference([full, sparse_pairs], labels, requested)
    @test isequal(average.distances, expected)
    @test average.counts == counts
    @test isnan(mean_dist([full, sparse_pairs]; sample_labels = labels,
        pairs = requested, nan_policy = :propagate).distances[1])
    @test isempty(mean_dist([sparse_pairs]; pairs = Tuple{Int, Int}[]).distances)
    @test mean_dist([sparse_pairs]; pairs = [(1, 2)]).counts == [1]
    missing_pair = (distances = Float64[], sample_labels = labels, pairs = Tuple{Int, Int}[])
    @test isnan(only(mean_dist([missing_pair]; pairs = [(1, 2)]).distances))
    # Missing pairs are unobserved and do not propagate a NaN.
    @test mean_dist([sparse_pairs, missing_pair]; sample_labels = labels,
        pairs = [(2, 3)], nan_policy = :propagate).distances == [8]
end

@testset "Mean distances: independent random averages and numerical range" begin
    rng = StableRNG(20261110)
    for case in 1:32
        n, draws = rand(rng, 2:10), rand(rng, 1:8)
        labels = ["s_$i" for i in randperm(rng, n)]
        inputs = NamedTuple[]
        for draw in 1:draws
            order = draw == 1 ? randperm(rng, n) : randperm(rng, n)[1:rand(rng, 0:n)]
            m = length(order)
            distances = zeros(m, m)
            for j in 1:m, i in j:m
                value = rand(rng) < 0.2 ? NaN : i == j ? 0.0 : rand(rng) * 10.0^rand(rng, -200:200)
                distances[i, j] = distances[j, i] = value
            end
            push!(inputs, (; distances, sample_labels = labels[order]))
        end
        requested = [(i, j) for j in 1:n for i in 1:n]
        for policy in (:skip, :propagate)
            expected, counts = mean_dist_reference(inputs, labels, requested; propagate = policy === :propagate)
            average = mean_dist(inputs; sample_labels = labels, nan_policy = policy)
            @test isequal(isnan.(vec(average.distances)), isnan.(expected))
            @test vec(average.counts) == counts
            valid = isfinite.(expected)
            @test all(isapprox.(vec(average.distances)[valid], expected[valid]; rtol = 8e-15, atol = 0))
            selected = vcat(requested, requested[1:2])
            result = mean_dist((input for input in inputs); sample_labels = labels,
                pairs = selected, nan_policy = policy)
            @test isequal(result.distances, vcat(vec(average.distances), vec(average.distances)[1:2]))
            @test result.counts == vcat(counts, counts[1:2])
        end
    end
    tiny, huge = nextfloat(0.0), floatmax(Float64)
    for values in ([huge, huge], [huge, 0.0, huge], [tiny, 0.0, tiny],
        [huge, 1.0, huge / 2], [0.0, 0.0, 0.0], [floatmin(Float64), tiny, 0.0])
        result = mean_dist(([0.0 x; x 0.0] for x in values))
        expected = Float64(sum(BigFloat.(values)) / length(values))
        @test isfinite(result.distances[1, 2])
        @test result.distances[1, 2] ≈ expected rtol=2e-15 atol=0
        @test result.counts[1, 2] == length(values)
    end
end

@testset "Mean distances: posterior beta pipeline and regionalization" begin
    trees = mean_dist_test_tree.([
        "((A:1,B:1):2,(C:2,D:2):1);",
        "((A:2,C:2):1,(B:1,D:1):2);",
        "((A:1,D:1):2,(B:2,C:2):1);",
        "((A:0,B:0):0,(C:0,D:0):0);",
    ])
    community = [1 1 0 0; 0 0 1 1; 1 0 1 0; 0 1 0 1]
    labels = ["west", "east", "middle", "edge"]
    requested = [(1, 2), (2, 1), (4, 3), (2, 2), (1, 2)]
    for family in (:sorensen, :jaccard), component in (:total, :turnover, :nestedness)
        inputs = [phylobeta(community, tree; sample_labels = labels, index_family = family) for tree in trees]
        average = mean_dist((input for input in inputs); component)
        plain = [(distances = getproperty(input, component), sample_labels = input.sample_labels) for input in inputs]
        expected, counts = mean_dist_reference(plain, labels, [(i, j) for j in 1:4 for i in 1:4])
        @test vec(average.distances) ≈ expected
        @test vec(average.counts) == counts
        @test average.component === component
        @test average.index_family === family
        @test average.replicates == 4
        paired = mean_dist((phylobeta(community, tree; sample_labels = labels,
            index_family = family, pairs = requested) for tree in trees); component, pairs = requested)
        @test paired.distances == [average.distances[i, j] for (i, j) in requested]
        @test paired.counts == [average.counts[i, j] for (i, j) in requested]
        direct = phyloregion(average; k = 2)
        baseline = phyloregion(average.distances; k = 2, sample_labels = labels)
        @test direct.membership == baseline.membership
        @test direct.sample_labels == labels
        @test direct.component === component
        context = prepare_phyloregion(average)
        @test context.component === component
        @test cutree(context; k = 2) == direct.membership
        @test select_linkage(average).correlations ≈ select_linkage(average.distances).correlations nans=true
        @test select_linkage(average).component === component
        optimal = optimal_phyloregion(average; k = 4)
        @test isequal(optimal.df, optimal_phyloregion(context; k = 4).df)
        @test optimal.component === component
        for call in (x -> phyloregion(x; k = 2), prepare_phyloregion, select_linkage, optimal_phyloregion)
            @test_throws ArgumentError call(paired)
        end
        poisoned = mean_dist(inputs; component, nan_policy = :propagate)
        @test_throws ArgumentError phyloregion(poisoned; k = 2)
    end
    taxonomic = mean_dist((beta_diss(community; sample_labels = labels) for _ in 1:3))
    @test taxonomic.distances == beta_diss(community).total
    @test all(==(3), taxonomic.counts)
    @test phyloregion(taxonomic; k = 2).sample_labels == labels
    empty_sites = mean_dist([beta_diss([1 0; 0 0])])
    @test empty_sites.distances[1, 2] == 1
    @test isnan(empty_sites.distances[2, 2])
    @test empty_sites.counts[2, 2] == 0
    @test_throws ArgumentError prepare_phyloregion(empty_sites)
end

@testset "Mean distances: validation, ownership and one-pass consumption" begin
    matrix = [0.0 2; 2 0]
    labels, requested = ["b", "a"], [(2, 1)]
    input = (distances = matrix, sample_labels = labels)
    result = mean_dist([input]; pairs = requested)
    matrix[1, 2] = matrix[2, 1] = 9
    labels[1] = "changed"
    requested[1] = (1, 1)
    @test result.distances == [2]
    @test result.sample_labels == ["b", "a"]
    @test result.pairs == [(2, 1)]
    result.distances[1] = 50
    @test matrix[1, 2] == 9
    draws = Ref(0)
    onepass = Iterators.Stateful((begin
        draws[] += 1
        [0.0 Float64(i); Float64(i) 0]
    end for i in 1:5))
    @test mean_dist(onepass).distances == [0 3; 3 0]
    @test draws[] == 5
    @test isempty(onepass)
    for bad in ([], [1], [ones(2, 3)], [ones(2, 2)], [[0.0 -1; -1 0]],
        [[0.0 Inf; Inf 0]], [[0.0 1; 2 0]], [[0.0 NaN; 0 0]],
        [(distances = [0 1; 1 0], sample_labels = ["a", "a"])],
        [(distances = [0 1; 1 0], sample_labels = ["a", " "])],
        [(distances = [0 1; 1 0], sample_labels = ["a"])],
        [(distances = [0 1; 1 0],)], [(sample_labels = ["a"],)],
        [(distances = [missing;;], sample_labels = ["a"])])
        @test_throws ArgumentError mean_dist(bad)
    end
    @test_throws ArgumentError mean_dist([matrix]; component = :other)
    @test_throws ArgumentError mean_dist([matrix]; nan_policy = :zero)
    @test_throws ArgumentError mean_dist([matrix]; sample_labels = "bad")
    @test_throws ArgumentError mean_dist([matrix]; sample_labels = ["a"])
    @test_throws ArgumentError mean_dist([matrix]; pairs = [(true, 1)])
    @test_throws ArgumentError mean_dist([matrix]; pairs = [(1, 3)])
    @test_throws ArgumentError mean_dist([input]; sample_labels = ["a", "b"])
    @test_throws ArgumentError mean_dist([zeros(2, 2), zeros(1, 1)])
    b1, b2 = beta_diss([1 0; 0 1]), beta_diss([1 0; 0 1]; index_family = :jaccard)
    @test_throws ArgumentError mean_dist([b1, b2])
    @test_throws ArgumentError mean_dist([b1, matrix])
    phylo = phylobeta([1 0; 0 1], mean_dist_test_tree("(a:1,b:1);"))
    @test_throws ArgumentError mean_dist([b1, phylo])
    for bad in (
        (distances = [1.0, 2], sample_labels = ["a", "b"], pairs = [(1, 2), (2, 1)]),
        (distances = [NaN, 2.0], sample_labels = ["a", "b"], pairs = [(1, 2), (2, 1)]),
        (distances = [1.0], sample_labels = ["a", "b"], pairs = [(1, 1)]),
        (distances = [Inf], sample_labels = ["a", "b"], pairs = [(1, 2)]),
        (distances = [-1.0], sample_labels = ["a", "b"], pairs = [(1, 2)]),
        (distances = [1.0], sample_labels = ["a", "b"], pairs = [(1, 3)]),
        (distances = Float64[], sample_labels = ["a", "b"], pairs = [(1, 2)]),
        (distances = zeros(2, 2), sample_labels = ["a", "b"], pairs = [(1, 2)]),
        (distances = [1.0], sample_labels = ["a", "b"]),
    )
        @test_throws ArgumentError mean_dist([bad]; pairs = [(1, 2)])
    end
    for value in (big"1e1000", big"1e-1000")
        @test_throws ArgumentError mean_dist([[big"0" value; value big"0"]])
    end
    @test_throws ArgumentError mean_dist([result]) # Do not give already averaged draws equal weight.
end

@testset "Mean distances: selected pairs use linear storage" begin
    labels = ["site_$i" for i in 1:30_000]
    selected = [(1, 30_000), (30_000, 1), (2, 2)]
    source = (distances = [0.75, 0.75, 0.0], sample_labels = labels, pairs = selected)
    mean_dist([source]; pairs = selected)
    allocated = @allocated result = mean_dist((source for _ in 1:2); pairs = selected)
    @test result.distances == [0.75, 0.75, 0]
    @test result.counts == [2, 2, 2]
    @test allocated < 80_000_000
    # Fail before allocating a 30,000-by-30,000 output without explicit pairs.
    @test_throws ArgumentError mean_dist([source])
    @test_throws ArgumentError mean_dist([(distances = [0.75], sample_labels = labels)])
end

function mean_dist_fixture_inputs(rows)
    inputs = NamedTuple[]
    for draw in groupby(rows, :draw)
        n = maximum(draw.i)
        distances = fill(NaN, n, n)
        labels = Vector{String}(undef, n)
        for row in eachrow(draw)
            distances[row.i, row.j] = distances[row.j, row.i] = row.distance
            labels[row.i], labels[row.j] = row.label_i, row.label_j
        end
        push!(inputs, (; distances, sample_labels = labels))
    end
    return inputs
end

@testset "Mean distances: saved local R mean_matrix and posterior beta references" begin
    fixture = CSV.read(joinpath(@__DIR__, "fixtures", "mean_dist", "r_reference.tsv"),
        DataFrame; delim = '\t')
    for rows in groupby(fixture, :case)
        labels = String.(split(first(rows.output_labels), ';'))
        inputs = mean_dist_fixture_inputs(rows)
        skip = mean_dist((input for input in inputs); sample_labels = labels)
        propagated = mean_dist(inputs; sample_labels = labels, nan_policy = :propagate)
        @test skip.replicates == length(inputs)
        @test skip.sample_labels == labels
        seen = Set{Tuple{Int, Int}}()
        requested = Tuple{Int, Int}[]
        for row in eachrow(rows)
            i, j = findfirst(==(row.label_i), labels), findfirst(==(row.label_j), labels)
            key = minmax(i, j)
            key in seen && continue
            push!(seen, key)
            push!(requested, (i, j))
            @test isapprox(skip.distances[i, j], row.mean_skip; nans = true, atol = 1e-13)
            @test skip.counts[i, j] == row.count
            @test isapprox(propagated.distances[i, j], row.mean_propagate; nans = true, atol = 1e-13)
            @test propagated.counts[i, j] == row.count
            if i != j && row.native_ok
                @test isapprox(propagated.distances[i, j], row.native_mean; nans = true, atol = 1e-13)
            end
        end
        # Convert sources to labeled pair vectors; repeated pairs still supply
        # only one observation per posterior draw.
        paired_inputs = map(inputs) do input
            n = length(input.sample_labels)
            pairs = [(i, j) for j in 1:n for i in j:n]
            push!(pairs, first(pairs))
            (distances = [input.distances[i, j] for (i, j) in pairs],
                sample_labels = input.sample_labels, pairs = pairs)
        end
        selected = mean_dist(paired_inputs; sample_labels = labels, pairs = requested)
        @test selected.distances ≈ [skip.distances[i, j] for (i, j) in requested] nans=true
        @test selected.counts == [skip.counts[i, j] for (i, j) in requested]
    end
end
