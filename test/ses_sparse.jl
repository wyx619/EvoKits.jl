include("ses_reference.jl")
const SES_MODELS = (:taxa_labels, :richness, :sample_pool, :phylogeny_pool,
    :frequency, :trialswap, :independentswap)

@testset "Sparse null models match dense draws and RNG consumption" begin
    fixtures = [rand(MersenneTwister(80), Bool, 7, 11), trues(4, 6), falses(4, 6),
        Bool[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0], trues(1, 4), trues(3, 1),
        falses(0, 4), falses(3, 0)]
    for presence in fixtures, model in SES_MODELS
        observed = [Int32.(findall(presence[i, :])) for i in axes(presence, 1)]
        saved = deepcopy(observed)
        workspace = EvoKits._CommunityNullWorkspace(observed, size(presence, 2))
        left, right = MersenneTwister(91), MersenneTwister(91)
        for _ in 1:5
            expected = ses_dense_draw(presence, model, 30, 50, right)
            rows = EvoKits._randomize_community!(workspace, observed, model, 30, 50, left)
            actual = fill(false, size(presence))
            for (site, row) in enumerate(rows)
                actual[site, row] .= true
                @test issorted(row) && allunique(row)
            end
            @test actual == expected
            @test rand(left, UInt64) == rand(right, UInt64)
            @test observed == saved
            @test sum(actual) == sum(presence)
            if model in (:richness, :sample_pool, :phylogeny_pool, :taxa_labels, :trialswap, :independentswap)
                @test sum(actual; dims = 2) == sum(presence; dims = 2)
            end
            if model in (:frequency, :trialswap, :independentswap)
                @test sum(actual; dims = 1) == sum(presence; dims = 1)
            elseif model === :taxa_labels
                @test actual * actual' == presence * presence'
            end
        end
    end
end

@testset "Sparse SES against independent dense metrics and null summaries" begin
    mktempdir() do dir
        path = joinpath(dir, "ses.nwk")
        # Integer branch lengths and weights make rank ties exactly comparable.
        write_tree(path, "((A:1,B:2):3,(C:4,D:1):2,E:3);")
        tree = read_tree(path)
        community = [2 0 3 0 1; 0 4 0 0 0; 0 0 0 0 0; 1 2 0 1 0; 1 1 1 1 1]
        for model in SES_MODELS
            for (metric, f, weighted, root) in ((:pd, ses_pd, false, true), (:pd, ses_pd, false, false),
                (:mpd, ses_mpd, false, true), (:mpd, ses_mpd, true, true),
                (:mntd, ses_mntd, false, true), (:mntd, ses_mntd, true, true))
                reference = ses_summary_reference(community, tree, metric, model; weighted, include_root = root)
                for input in (community, sparse(community))
                    options = metric === :pd ? (; include_root = root) : (; abundance_weighted = weighted)
                    result = f(input, tree; null_model = model, runs = 11, iterations = 20,
                        max_attempts_per_swap = 40, rng = MersenneTwister(123), options...)
                    prefix = String(metric)
                    @test same_values(getproperty(result, Symbol(prefix * "_obs")), reference.obs)
                    @test same_values(getproperty(result, Symbol(prefix * "_rand_mean")), reference.means)
                    @test same_values(getproperty(result, Symbol(prefix * "_rand_sd")), reference.sds)
                    @test same_values(getproperty(result, Symbol(prefix * "_obs_rank")), reference.ranks)
                    @test same_values(getproperty(result, Symbol(prefix * "_obs_p")), reference.p)
                    actual_z = getproperty(result, Symbol(prefix * "_obs_z"))
                    # Small rounding in online variance is immaterial unless
                    # the reference SD is exactly zero, when both must be NaN.
                    @test same_values(actual_z, reference.z)
                    @test result.ntaxa == vec(sum(community .> 0; dims = 2))
                    @test result.null_model === model
                    @test result.runs == 11
                    if metric !== :pd
                        @test same_values(getproperty(result, metric === :mpd ? :mpd_nri : :mntd_nti), -actual_z)
                    end
                end
            end
        end
        # Repeating a constant community must retain PD in every draw. The old
        # per-site marker reused between draws incorrectly skipped its branches.
        result = ses_pd(ones(1, 5), tree; runs = 4, rng = MersenneTwister(1))
        @test result.pd_obs == result.pd_rand_mean == [sum(tree.edge_length)]
        @test result.pd_rand_sd == [0]
        @test isnan(result.pd_obs_z[1])
        @test result.pd_obs_rank == [3]
        @test result.pd_obs_p == [3 / 5]
        for f in (ses_pd, ses_mpd, ses_mntd)
            for options in ((runs = 0,), (iterations = -1,), (max_attempts_per_swap = 0,), (null_model = :wrong,))
                @test_throws ArgumentError f(community, tree; options...)
            end
            @test_throws ArgumentError f(sparse([-1 0 0 0 0]), tree)
            @test_throws ArgumentError f(community, tree; sample_labels = ["bad"])
            @test_throws ArgumentError f(community, tree; tip_labels = ["A", "B", "C", "D", "X"])
            @test isempty(getproperty(f(spzeros(0, 5), tree; runs = 2),
                f === ses_pd ? :pd_obs : f === ses_mpd ? :mpd_obs : :mntd_obs))
        end
        # No finite null values are reported as NaN rather than an artificial 0.
        @test isnan(ses_pd([1 0 0 0 0], tree; runs = 3, include_root = false).pd_rand_mean[1])
        @test isnan(ses_mpd([1 0 0 0 0], tree; runs = 3).mpd_rand_mean[1])
        @test isnan(ses_mntd([1 0 0 0 0], tree; runs = 3).mntd_rand_mean[1])
    end
end

@testset "Sparse SES input order, snapshot and zero-root cancellation" begin
    rng = MersenneTwister(20261109)
    mktempdir() do dir
        path = joinpath(dir, "random.nwk")
        for trial in 1:12
            p, n = rand(rng, 3:20), rand(rng, 2:8)
            write_tree(path, random_community_tree(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:4, n, p)
            matrix[rand(rng, n, p) .< 0.7] .= 0
            permutation = randperm(rng, p)
            reordered = sparse(matrix[:, permutation])
            saved = copy(reordered)
            labels = ["site_$i" for i in 1:n]
            data = CommunityData(reordered; sample_labels = labels, species_labels = tree.tip_labels[permutation])
            for (f, weighted) in ((ses_pd, false), (ses_mpd, false), (ses_mpd, true), (ses_mntd, false), (ses_mntd, true))
                options = f === ses_pd ? (;) : (; abundance_weighted = weighted)
                a = f(matrix, tree; runs = 7, rng = MersenneTwister(trial), sample_labels = labels, options...)
                b = f(reordered, tree; runs = 7, rng = MersenneTwister(trial),
                    sample_labels = labels, tip_labels = tree.tip_labels[permutation], options...)
                c = f(data, tree; runs = 7, rng = MersenneTwister(trial), options...)
                check_context_result(a, b)
                check_context_result(a, c)
            end
            @test reordered == saved && nonzeros(reordered) == nonzeros(saved)
        end
        write_tree(path, "((A:0.000001,B:0.000002):1000000000000,C:1);")
        tree = read_tree(path)
        result = ses_pd([1 1 0], tree; null_model = :trialswap, iterations = 0,
            runs = 3, include_root = false)
        @test result.pd_obs[1] ≈ 3e-6
        @test result.pd_rand_mean[1] ≈ 3e-6
        @test result.pd_rand_sd[1] == 0
        bad = deepcopy(tree)
        bad.edge_length[1] = -1
        for f in (ses_pd, ses_mpd, ses_mntd)
            @test_throws ArgumentError f([1 1 0], bad)
        end
        # Virtual edges may span multiple original branches. Their lengths
        # must survive cancellation beneath a much longer common root path.
        for (newick, labels, expected_distance) in (
            ("((A:0.000001,B:0.000002):1000000000000,C:1);", ["A", "B", "C"], 3e-6),
            ("(((A:0.000001,D:1):0.000002,B:0.000003):1000000000000,C:1);", ["A", "B", "C", "D"], 6e-6))
            write_tree(path, newick)
            precision_tree = read_tree(path)
            community = reshape([1, 2, zeros(Int, length(labels) - 2)...], 1, :)
            for (metric, f) in ((:mpd, ses_mpd), (:mntd, ses_mntd)), weighted in (false, true)
                result = f(sparse(community), precision_tree; tip_labels = labels,
                    null_model = :trialswap, iterations = 0, runs = 3, abundance_weighted = weighted)
                expected = metric === :mpd && weighted ? 4expected_distance / 9 : expected_distance
                @test getproperty(result, Symbol(metric, :_obs))[1] ≈ expected
                @test getproperty(result, Symbol(metric, :_rand_mean))[1] ≈ expected
                @test getproperty(result, Symbol(metric, :_rand_sd))[1] == 0
            end
        end
    end
end

@testset "SES sparse workspace memory has no site-by-species buffers" begin
    rng = MersenneTwister(20261110)
    n, p = 2000, 8000
    community = sparse(rand(rng, 1:n, 6000), rand(rng, 1:p, 6000), ones(Int, 6000), n, p)
    observed = EvoKits._community_tip_rows(community, 1:p)
    workspace = EvoKits._CommunityNullWorkspace(observed, p)
    @test Base.summarysize(workspace) < 2_000_000
    for model in SES_MODELS
        EvoKits._randomize_community!(workspace, observed, model, 10, 20, rng)
        @test Base.summarysize(workspace) < 2_000_000
    end
    mktempdir() do dir
        path = joinpath(dir, "large.nwk")
        write_tree(path, "(" * join(["s$i:1" for i in 1:p], ',') * ");")
        tree = read_tree(path)
        for (f, weighted) in ((ses_pd, false), (ses_mpd, false), (ses_mpd, true), (ses_mntd, true))
            options = f === ses_pd ? (;) : (; abundance_weighted = weighted)
            f(community, tree; runs = 1, options...)
            bytes = @allocated f(community, tree; runs = 1, options...)
            # One old Bool matrix was 16 MB, before random presence, temporary
            # buffers and (for weighted metrics) a 128 MB abundance matrix.
            @test bytes < 12_000_000
        end
    end
end
