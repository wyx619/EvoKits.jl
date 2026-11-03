using Test, EvoKits, Random, StableRNGs, Statistics, SparseArrays, DataFrames, CSV
include("beta_ses_reference.jl")

const BETA_SES_MODELS = (:taxa_labels, :richness, :sample_pool, :phylogeny_pool,
    :frequency, :trialswap, :independentswap)
const BETA_SES_ARRAYS = (:observed, :null_mean, :null_sd, :z_score, :p_value, :valid_replicates)

function beta_ses_test_tree(newick)
    return mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, newick)
        read_tree(path)
    end
end

@testset "Beta SES: independent dense null draws, BigFloat metrics and finite counts" begin
    tree = beta_ses_test_tree("((A:1,B:2):3,(C:4,D:1):2,E:0);")
    community = [2 0 3 0 1; 0 4 0 0 0; 0 0 0 0 0; 1 2 0 1 0; 1 1 1 1 1]
    requested = [(1, 2), (2, 1), (3, 4), (3, 3), (5, 2), (1, 2), (4, 4)]
    for model in BETA_SES_MODELS, family in (:sorensen, :jaccard), component in (:total, :turnover, :nestedness)
        reference = beta_ses_reference(community, tree, model, family, component)
        full_rng, pair_rng = StableRNG(42), StableRNG(42)
        full = phylobeta_ses(community, tree; null_model = model, index_family = family,
            component, runs = 13, iterations = 15, max_attempts_per_swap = 35, rng = full_rng)
        selected = phylobeta_ses(sparse(community), tree; null_model = model, index_family = family,
            component, runs = 13, iterations = 15, max_attempts_per_swap = 35, rng = pair_rng, pairs = requested)
        for field in BETA_SES_ARRAYS
            actual, expected = getproperty(full, field), getproperty(reference, field)
            if field in (:p_value, :valid_replicates)
                @test isequal(vec(actual), expected)
            elseif field === :z_score
                # Both references agree on constant/null boundaries; random
                # means/SDs may differ in their last bits due to online order.
                @test isapprox(vec(actual), expected; atol = 1e-10, rtol = 1e-9, nans = true)
            else
                @test isapprox(vec(actual), expected; atol = 1e-13, rtol = 1e-12, nans = true)
            end
            @test isequal(getproperty(selected, field), [actual[i, j] for (i, j) in requested])
            @test isequal(actual, transpose(actual))
        end
        @test rand(full_rng, UInt64) == reference.next_random
        @test rand(pair_rng, UInt64) == reference.next_random
        @test full.richness == [3, 1, 0, 3, 5]
        @test full.null_model === model
        @test full.index_family === family && full.component === component
        @test selected.pairs == requested
        @test full.runs == 13 && full.iterations == 15 && full.max_attempts_per_swap == 35
    end
end

@testset "Beta SES: no valid draws, one draw, constant nulls and empty observations" begin
    tree = beta_ses_test_tree("((a:1,b:1):1,(c:2,d:2):0);")
    same = [1 0 0 0; 0 1 0 0; 0 0 0 0]
    for component in (:total, :turnover, :nestedness)
        constant = phylobeta_ses(same, tree; component, null_model = :trialswap, iterations = 0, runs = 4)
        defined = isfinite.(constant.observed)
        @test isequal(constant.observed, constant.null_mean)
        @test all(==(0), constant.null_sd[defined])
        @test all(isnan, constant.z_score)
        @test all(==(4), constant.valid_replicates[defined])
        @test all(==(3 / 5), constant.p_value[defined])
        @test all(isnan, constant.p_value[.!defined])
        one = phylobeta_ses(same, tree; component, null_model = :trialswap, iterations = 0, runs = 1)
        @test all(isnan, one.null_sd)
        @test all(isnan, one.z_score)
        @test all(==(0.75), one.p_value[defined])
        zero_tree = beta_ses_test_tree("(a:0,b:0,c:0,d:0);")
        undefined = phylobeta_ses(same, zero_tree; component, runs = 3, rng = StableRNG(12))
        for field in BETA_SES_ARRAYS[1:5]
            @test all(isnan, getproperty(undefined, field))
        end
        @test all(iszero, undefined.valid_replicates)
    end
    # Frequency randomization can move occurrences into originally empty sites;
    # finite null statistics remain reportable when observed beta is undefined.
    result = phylobeta_ses([1 1 0 0; 0 0 1 1; 0 0 0 0], tree;
        null_model = :frequency, component = :turnover, runs = 30, rng = StableRNG(45))
    @test isnan(result.observed[1, 3])
    @test result.valid_replicates[1, 3] > 0
    @test isfinite(result.null_mean[1, 3])
    @test isnan(result.z_score[1, 3]) && isnan(result.p_value[1, 3])
    @test 0 < result.valid_replicates[1, 2] < 30
    @test isnan(result.observed[3, 3]) && isnan(result.p_value[3, 3])
end

@testset "Beta SES: fixed local R beta and online-summary references" begin
    fixture = CSV.read(joinpath(@__DIR__, "fixtures", "beta_ses", "r_reference.tsv"), DataFrame; delim = '\t')
    for rows in groupby(fixture, [:case, :model, :family, :component])
        config = first(eachrow(rows))
        tree = beta_ses_test_tree(config.newick)
        community = reduce(vcat, [permutedims(parse.(Int, split(row, ','))) for row in split(config.community, ';')])
        pairs = [(Int(row.i), Int(row.j)) for row in eachrow(rows)]
        result = phylobeta_ses(sparse(community), tree; null_model = Symbol(config.model),
            index_family = Symbol(config.family), component = Symbol(config.component), runs = config.runs,
            iterations = config.iterations, max_attempts_per_swap = config.attempts, pairs, rng = StableRNG(config.seed))
        for field in (:observed, :null_mean, :null_sd, :z_score)
            @test isapprox(getproperty(result, field), rows[!, field]; nans = true, atol = 1e-11, rtol = 1e-9)
        end
        @test result.valid_replicates == rows.valid_replicates
        @test result.pairs == pairs
        @test result.runs == config.runs
        @test all(eachindex(result.p_value)) do i
            isnan(result.p_value[i]) ? isnan(rows.p_min[i]) : rows.p_min[i] - 1e-14 <= result.p_value[i] <= rows.p_max[i] + 1e-14
        end
    end
end

@testset "Beta SES: label alignment, full tree pool, aliases and ownership" begin
    tree = beta_ses_test_tree("((a:1,b:2):3,(c:4,d:1):2,e:3);")
    matrix = [1 0 2 0 0; 0 3 0 0 0; 0 0 0 0 0; 1 1 0 0 0]
    labels, order = ["west", "east", "empty", "mixed"], [3, 1, 5, 2, 4]
    data = CommunityData(sparse(matrix[:, order]); sample_labels = labels, species_labels = tree.tip_labels[order])
    subset = CommunityData(sparse(matrix[:, [1, 2, 3]]); sample_labels = labels, species_labels = ["a", "b", "c"])
    for model in BETA_SES_MODELS
        a = phylobeta_ses(matrix, tree; sample_labels = labels, null_model = model, runs = 9, rng = StableRNG(8))
        b = phylobeta_ses(sparse(matrix[:, order]), tree; sample_labels = labels,
            tip_labels = tree.tip_labels[order], null_model = model, runs = 9, rng = StableRNG(8))
        c = phylobeta_ses(data, tree; null_model = model, runs = 9, rng = StableRNG(8))
        d = phylobeta_ses(subset, tree; null_model = model, runs = 9, rng = StableRNG(8))
        for field in BETA_SES_ARRAYS
            @test isequal(getproperty(a, field), getproperty(b, field))
            @test isequal(getproperty(a, field), getproperty(c, field))
            @test isequal(getproperty(a, field), getproperty(d, field))
        end
        @test a.sample_labels == labels
    end
    a = phylobeta_ses(data, tree; runs = 9, rng = StableRNG(8))
    for alias in (:tipshuffle, Symbol("taxa.labels"))
        b = phylobeta_ses(data, tree; null_model = alias, runs = 9, rng = StableRNG(8))
        @test isequal(a.null_mean, b.null_mean)
        @test b.null_model === :taxa_labels
    end
    saved_tree, saved_matrix = deepcopy(tree), copy(data.matrix)
    requested = [(1, 2), (2, 1)]
    result = phylobeta_ses(data, tree; pairs = requested, runs = 3, rng = StableRNG(9))
    requested[1] = (1, 1)
    labels[1] = "changed"
    result.observed[1] = 50
    result.sample_labels[1] = "modified"
    @test result.pairs == [(1, 2), (2, 1)]
    @test data.matrix == saved_matrix
    @test data.sample_labels[1] == "west"
    @test tree.edge_length == saved_tree.edge_length && tree.tip_labels == saved_tree.tip_labels
    @test phylobeta_ses(data, tree; pairs = [(1, 2)], runs = 3, rng = StableRNG(9)).observed[1] != 50
    # Reordering sites preserves non-random observed values. Random richness
    # and swap models depend on input order; no invariance is promised for them.
    reordered = phylobeta_ses(matrix[[4, 1, 2, 3], :], tree; runs = 9, rng = StableRNG(8))
    @test isequal(reordered.observed, a.observed[[4, 1, 2, 3], [4, 1, 2, 3]])
    @test isequal(reordered.null_mean, a.null_mean[[4, 1, 2, 3], [4, 1, 2, 3]])
end

@testset "Beta SES: sparse scale, extreme branch lengths and invalid inputs" begin
    tree = beta_ses_test_tree("((a:1,b:2):3,(c:4,d:1):2,e:3);")
    matrix = [1 0 2 0 0; 0 3 0 0 0; 0 0 0 0 0; 1 1 0 0 0]
    for component in (:total, :turnover, :nestedness), scale in (1e-200, 1e200)
        scaled = deepcopy(tree)
        scaled.edge_length .*= scale
        result = phylobeta_ses(matrix, scaled; component, runs = 7, rng = StableRNG(3))
        baseline = phylobeta_ses(matrix, tree; component, runs = 7, rng = StableRNG(3))
        # Distances and standardized scores are invariant to tree-wide scale;
        # exact rank ties are defined on computed Float64 values, not a tolerance.
        for field in (:observed, :null_mean, :null_sd, :z_score)
            @test isapprox(getproperty(result, field), getproperty(baseline, field);
                nans = true, atol = 1e-12, rtol = 1e-10)
        end
        @test result.valid_replicates == baseline.valid_replicates
    end
    long = beta_ses_test_tree("(((a:1,b:1):1):1e20,c:1);")
    @test phylobeta_ses([1 0 0; 0 1 0], long; null_model = :trialswap, iterations = 0,
        runs = 2, pairs = [(1, 2)]).observed[1] > 0
    for options in ((runs = 0,), (runs = true,), (runs = 1.5,), (iterations = -1,),
        (iterations = false,), (max_attempts_per_swap = 0,), (max_attempts_per_swap = 1.5,),
        (null_model = :wrong,), (null_model = :rowwise,), (null_model = :colwise,),
        (index_family = :simpson,), (component = :wrong,), (pairs = [(1, 6)],),
        (pairs = [(true, 1)],), (sample_labels = ["one"],),
        (sample_labels = ["x", "x", "y", "z"],), (sample_labels = [" ", "x", "y", "z"],),
        (tip_labels = ["a", "a", "c", "d", "e"],), (tip_labels = ["a", "b", "c", "d", "X"],))
        @test_throws ArgumentError phylobeta_ses(matrix, tree; runs = 2, options...)
    end
    for bad in ([-1 0 0 0 0], [NaN 0 0 0 0], sparse([Inf 0 0 0 0]), zeros(2, 3))
        @test_throws ArgumentError phylobeta_ses(bad, tree; runs = 2)
    end
    for value in (-1.0, Inf, NaN)
        bad = deepcopy(tree)
        bad.edge_length[1] = value
        @test_throws ArgumentError phylobeta_ses(matrix, bad; runs = 2)
    end
    data = CommunityData(matrix; species_labels = tree.tip_labels)
    @test_throws ArgumentError phylobeta_ses(data, tree; sample_labels = data.sample_labels)
    @test_throws ArgumentError phylobeta_ses(data, tree; tip_labels = tree.tip_labels)
    unknown = CommunityData(zeros(2, 1); species_labels = ["unknown"])
    @test_throws ArgumentError phylobeta_ses(unknown, tree; runs = 2)
    for model in BETA_SES_MODELS
        empty_result = phylobeta_ses(spzeros(0, 5), tree; runs = 2, null_model = model)
        @test size(empty_result.observed) == (0, 0)
        @test isempty(empty_result.richness)
        no_pairs = phylobeta_ses(matrix, tree; pairs = [], null_model = model, runs = 2)
        @test size(no_pairs.observed) == (0,)
    end
    @test phylobeta_ses(matrix, tree; runs = 2, sample_labels = [1, 2, 3, 4]).sample_labels == ["1", "2", "3", "4"]
end

@testset "Beta SES: memory scales with occurrences and requested pairs" begin
    n, p = 4000, 10000
    rng = StableRNG(71)
    matrix = sparse(rand(rng, 1:n, 12000), rand(rng, 1:p, 12000), ones(Int, 12000), n, p)
    tree = beta_ses_test_tree("(" * join(["sp_$i:1" for i in 1:p], ',') * ");")
    requested = [(1, n), (n, 1), (2, 2), (1, n)]
    phylobeta_ses(matrix, tree; runs = 2, pairs = requested, rng = StableRNG(7))
    bytes = @allocated result = phylobeta_ses(matrix, tree; runs = 3, pairs = requested, rng = StableRNG(7))
    @test bytes < 15_000_000 # One dense Bool community alone would need 40 MB.
    @test result.observed[1] == result.observed[2] == result.observed[4]
    @test Base.summarysize(result) < 1_000_000
    # Many replicates retain only online arrays. No pairs-by-runs history.
    small_tree = beta_ses_test_tree("(a:1,b:1,c:1,d:1);")
    small = [1 1 0 0; 0 1 1 0; 0 0 0 1]
    many_pairs = repeat([(1, 2), (2, 3), (1, 1)], 1000)
    phylobeta_ses(small, small_tree; runs = 2, pairs = many_pairs, rng = StableRNG(7))
    bytes = @allocated phylobeta_ses(small, small_tree; runs = 300, pairs = many_pairs, rng = StableRNG(7))
    @test bytes < 3_000_000 # Saving the 900,000 draws would alone take 7.2 MB.
end
