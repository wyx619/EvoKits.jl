using Statistics

const REGION_TEST_METHODS = (:average, :single, :complete, :ward_D, :ward_D2,
    :mcquitty, :median, :centroid)

# Independent full active-distance dictionary in BigFloat. At each emitted
# merge verify the global minimum, then update by the defining recurrence.
# No production packed indexing, nearest-neighbor cache or chain is used.
function check_linkage_recurrence(distances, result)
    n = size(distances, 1)
    active = Dict(-i => [i] for i in 1:n)
    current = Dict(minmax(-i, -j) => BigFloat(distances[i, j]) for j in 2:n for i in 1:(j - 1))
    method = result.method
    method === :ward_D2 && (current = Dict(pair => value^2 for (pair, value) in current))
    for step in 1:(n - 1)
        a, b = result.merge[step, 1], result.merge[step, 2]
        @test a != b && a < step && b < step && haskey(active, a) && haskey(active, b)
        d = current[minmax(a, b)]
        best = minimum(values(current))
        expected = method === :ward_D2 ? sqrt(d) : d
        @test result.heights[step] ≈ Float64(expected) rtol = 1e-11 atol = 1e-13
        @test Float64(d) ≈ Float64(best) rtol = 1e-11 atol = 1e-13
        left_size, right_size = length(active[a]), length(active[b])
        updates = Dict{Tuple{Int, Int}, BigFloat}()
        for other in keys(active)
            other in (a, b) && continue
            x, y = current[minmax(a, other)], current[minmax(b, other)]
            z = length(active[other])
            p, q = BigFloat(left_size), BigFloat(right_size)
            value = if method === :single
                min(x, y)
            elseif method === :complete
                max(x, y)
            elseif method === :average
                (p * x + q * y) / (p + q)
            elseif method === :mcquitty
                (x + y) / 2
            elseif method in (:ward_D, :ward_D2)
                ((p + z) * x + (q + z) * y - z * d) / (p + q + z)
            elseif method === :median
                (x + y) / 2 - d / 4
            else
                (p * x + q * y) / (p + q) - p * q * d / (p + q)^2
            end
            updates[minmax(step, other)] = value
        end
        filter!(pair -> !(a in pair.first || b in pair.first), current)
        merge!(current, updates)
        active[step] = vcat(active[a], active[b])
        delete!(active, a)
        delete!(active, b)
    end
    @test length(active) == 1 && sort!(only(values(active))) == collect(1:n)
end

function linkage_reference_cophenetic(result)
    n = length(result.sample_labels)
    descendants = Dict(-i => [i] for i in 1:n)
    values = zeros(n, n)
    for step in axes(result.merge, 1)
        a, b = result.merge[step, 1], result.merge[step, 2]
        for i in descendants[a], j in descendants[b]
            values[i, j] = values[j, i] = result.heights[step]
        end
        descendants[step] = vcat(descendants[a], descendants[b])
    end
    return values
end

function linkage_reference_ss(distances, membership)
    n = length(membership)
    tss = sum(BigFloat(distances[i, j])^2 for j in 2:n for i in 1:(j - 1); init = BigFloat(0)) / n
    wss = BigFloat(0)
    for region in unique(membership)
        sites = findall(==(region), membership)
        wss += sum(BigFloat(distances[sites[i], sites[j]])^2 for j in 2:length(sites)
            for i in 1:(j - 1); init = BigFloat(0)) / length(sites)
    end
    return Float64(tss), Float64(wss), tss == 0 ? NaN : Float64(1 - wss / tss)
end

@testset "Eight linkage methods: independent global recurrences and inversions" begin
    inversion = [0 1 1.1; 1 0 1.2; 1.1 1.2 0]
    for method in (:median, :centroid)
        result = phyloregion(inversion; k = 2, method)
        @test result.heights ≈ [1, 0.9]
        @test result.merge == [-1 -2; 1 -3]
        @test cutree(result; k = 2) == [1, 1, 2]
    end
    rng = MersenneTwister(2041)
    for trial in 1:20
        n = rand(rng, 3:12)
        raw = trial <= 10 ? Float64.(rand(rng, 0:6, n, n)) : rand(rng, n, n)
        distances = Matrix(Symmetric(raw))
        distances[diagind(distances)] .= 0
        for method in REGION_TEST_METHODS
            result = phyloregion(distances; k = rand(rng, 1:n), method)
            @test result.method === method
            check_linkage_recurrence(distances, result)
            expected = region_reference_cut(result.merge, n, result.k)
            regions, ed = region_reference_summary(distances, expected)
            @test result.membership == expected
            @test same_values(result.region_distances, regions)
            @test same_values(result.distinctiveness, ed)
        end
    end
    for method in REGION_TEST_METHODS
        @test isempty(phyloregion(zeros(1, 1); k = 1, method).heights)
        for factor in (1e-300, 1e300)
            result = phyloregion(inversion .* factor; k = 2, method)
            baseline = phyloregion(inversion; k = 2, method)
            @test result.membership == baseline.membership
            @test result.heights ./ factor ≈ baseline.heights
        end
    end
    @test phyloregion(inversion; k = 1, method = Symbol("ward.D")).method === :ward_D
    @test phyloregion(inversion; k = 1, method = Symbol("ward.D2")).method === :ward_D2
    @test_throws ArgumentError phyloregion([0 0 1.7e308 1.7e308; 0 0 1.7e308 1.7e308;
        1.7e308 1.7e308 0 0; 1.7e308 1.7e308 0 0];
        k = 1, method = :ward_D)
end

@testset "Prepared regionalization: reusable cuts, beta inputs and ownership" begin
    distances = [0 1 8 10; 1 0 6 12; 8 6 0 2; 10 12 2 0]
    labels = ["a", "b", "c", "d"]
    for method in REGION_TEST_METHODS
        context = prepare_phyloregion(distances; sample_labels = labels, method)
        @test context isa PhyloregionContext
        @test context.sample_labels == labels && context.method === method
        @test length(context.distances) == 6
        for k in 1:4
            result = phyloregion(context; k)
            direct = phyloregion(distances; k, sample_labels = labels, method)
            @test result.membership == direct.membership == cutree(context; k) == cutree(direct; k)
            @test result.merge == direct.merge && result.heights == direct.heights
            @test result.region_distances == direct.region_distances
            @test same_values(result.distinctiveness, direct.distinctiveness)
            @test result.sample_labels == labels && result.component === :distance
        end
        snapshot = deepcopy(context)
        result = phyloregion(context; k = 2)
        result.membership .= 99
        result.merge .= 99
        result.heights .= 99
        result.sample_labels[1] = "modified"
        result.region_distances .= 99
        result.cluster_sizes .= 99
        result.distinctiveness .= 99
        cutree(context; k = 2) .= 99
        @test context.merge == snapshot.merge && context.heights == snapshot.heights
        @test context.sample_labels == snapshot.sample_labels && context.distances == snapshot.distances
        for invalid in (0, 5, true, 1.5)
            @test_throws ArgumentError cutree(context; k = invalid)
            @test_throws ArgumentError phyloregion(context; k = invalid)
        end
    end
    context = prepare_phyloregion(distances; sample_labels = labels)
    distances .= 0
    labels[1] = "modified"
    @test context.heights == [1, 2, 9] && context.sample_labels[1] == "a"
    @test phyloregion(context; k = 2).region_distances == [0 9; 9 0]
    data = CommunityData(sparse([1 1 0; 1 0 1; 0 1 1]); sample_labels = ["x", "y", "z"])
    beta = beta_diss(data)
    for component in (:total, :turnover, :nestedness)
        prepared = prepare_phyloregion(beta; component)
        @test prepared.component === component && prepared.sample_labels == data.sample_labels
        @test cutree(prepared; k = 2) == phyloregion(beta; k = 2, component).membership
        @test select_linkage(prepared).component === component
        @test optimal_phyloregion(prepared).component === component
    end
    @test_throws ArgumentError prepare_phyloregion(beta; component = :wrong)
    @test_throws ArgumentError prepare_phyloregion(beta_diss(data; pairs = [(1, 2)]))
    @test_throws ArgumentError prepare_phyloregion([0 NaN; NaN 0])
end

@testset "Linkage correlations and k profiles: independent original-pair formulas" begin
    rng = MersenneTwister(665)
    for trial in 1:12
        n = rand(rng, 3:12)
        distances = Matrix(Symmetric(rand(rng, n, n)))
        distances[diagind(distances)] .= 0
        selection = select_linkage(distances)
        @test selection isa LinkageSelectionResult
        @test selection.methods == collect(REGION_TEST_METHODS)
        @test selection.best_method === selection.methods[argmax(selection.correlations)]
        pairs = [(i, j) for j in 2:n for i in 1:(j - 1)]
        for (index, method) in enumerate(REGION_TEST_METHODS)
            context = prepare_phyloregion(distances; method)
            result = phyloregion(context; k = 2)
            cophenetic = linkage_reference_cophenetic(result)
            expected = cor([distances[i, j] for (i, j) in pairs], [cophenetic[i, j] for (i, j) in pairs])
            @test selection.correlations[index] ≈ expected
            @test same_values(select_linkage(context).correlations, selection.correlations)
            profile = optimal_phyloregion(context; k = n)
            @test profile isa OptimalPhyloregionResult
            @test profile.method === method && profile.component === :distance
            @test profile.df.k == collect(1:n)
            for k in 1:n
                membership = region_reference_cut(context.merge, n, k)
                total, within, ev = linkage_reference_ss(distances, membership)
                @test profile.df.tss[k] ≈ total
                @test profile.df.totwss[k] ≈ within atol = 1e-13
                @test profile.df.totbss[k] ≈ total - within atol = 1e-13
                @test profile.df.ev[k] ≈ ev atol = 1e-13
            end
            @test profile.df.ev[1] == 0 && profile.df.ev[end] == 1
            @test optimal_phyloregion(context; k = 2).df.ev == profile.df.ev[1:2]
        end
    end
end

@testset "Regionalization selection: undefined scores, thresholds and scaling" begin
    distances = [0 1 8 10; 1 0 6 12; 8 6 0 2; 10 12 2 0]
    selection = select_linkage(distances; methods = (:single, :average))
    @test selection.methods == [:single, :average]
    for matrix in (zeros(1, 1), zeros(4, 4), ones(4, 4) - Matrix{Float64}(I, 4, 4))
        @test all(isnan, select_linkage(matrix).correlations)
        @test select_linkage(matrix).best_method === nothing
    end
    # With only one original pair Pearson correlation is undefined.
    @test select_linkage([0 1; 1 0]).best_method === nothing
    for matrix in (zeros(1, 1), zeros(4, 4))
        result = optimal_phyloregion(matrix)
        @test result.optimal_k === nothing && result.reason === :zero_dispersion
        @test all(isnan, result.df.ev) && all(iszero, result.df.tss)
    end
    @test optimal_phyloregion(distances; k = 1).reason === :thresholds_not_met
    @test optimal_phyloregion(distances; k = 4, inc_thresholds = [0.1], ev_thresholds = [0.9]).optimal_k == 2
    @test optimal_phyloregion(distances; k = 4, inc_thresholds = [1], ev_thresholds = [1]).optimal_k == 4
    @test optimal_phyloregion(distances; k = 4, inc_thresholds = [0], ev_thresholds = [1], precision = 15).optimal_k === nothing
    for methods in ((), (:average, :average), (:ward_D, Symbol("ward.D")), (:wrong,), ("average",), :average)
        @test_throws ArgumentError select_linkage(distances; methods)
    end
    for thresholds in ([], [-0.1], [1.1], [NaN], [Inf], [true])
        @test_throws ArgumentError optimal_phyloregion(distances; inc_thresholds = thresholds)
        @test_throws ArgumentError optimal_phyloregion(distances; ev_thresholds = thresholds)
    end
    for precision in (-1, 16, true, 1.5)
        @test_throws ArgumentError optimal_phyloregion(distances; precision)
    end
    for k in (0, 5, true, 1.5)
        @test_throws ArgumentError optimal_phyloregion(distances; k)
    end
    @test optimal_phyloregion(distances; inc_thresholds = [0.1, 0.01, 0.1], ev_thresholds = [0.5, 0.9]).optimal_k == 3
    for factor in (1e-300, 1e300)
        scaled = Float64.(distances) .* factor
        @test same_values(select_linkage(scaled).correlations, select_linkage(distances).correlations)
        result = optimal_phyloregion(scaled; k = 4)
        baseline = optimal_phyloregion(distances; k = 4)
        @test same_values(result.df.ev, baseline.df.ev)
        @test result.optimal_k == baseline.optimal_k
    end
end

@testset "Eight linkages and k selection: saved local R outputs" begin
    readfixture(name) = CSV.read(joinpath(@__DIR__, "fixtures", name), DataFrame;
        delim = '\t', missingstring = "NaN")
    inputs = readfixture("regionalization_linkage_input.tsv")
    merges = readfixture("regionalization_linkage_merges_r.tsv")
    cuts = readfixture("regionalization_linkage_cuts_r.tsv")
    scores = readfixture("regionalization_linkage_scores_r.tsv")
    profiles = readfixture("regionalization_optimal_r.tsv")
    for case in unique(inputs.case)
        input = inputs[inputs.case .== case, :]
        n = maximum(input.i)
        distances = reshape(input.distance, n, n)
        selection = select_linkage(distances)
        for (index, method) in enumerate(REGION_TEST_METHODS)
            r_method = method === :ward_D ? "ward.D" : method === :ward_D2 ? "ward.D2" : string(method)
            r_merge = merges[(merges.case .== case) .& (merges.method .== r_method), :]
            r_score = scores[(scores.case .== case) .& (scores.method .== r_method), :]
            r_profile = profiles[(profiles.case .== case) .& (profiles.method .== r_method), :]
            context = prepare_phyloregion(distances; method)
            @test same_values(context.heights, r_merge.height)
            @test selection.correlations[index] ≈ only(r_score.correlation)
            result = optimal_phyloregion(context; k = n)
            @test same_values(result.df.ev, r_profile.ev)
            @test same_values(result.df.totwss, r_profile.totwss)
            @test same_values(result.df.totbss, r_profile.totbss)
            @test same_values(result.df.tss, r_profile.tss)
            if ismissing(r_profile.optimal_k[1])
                @test result.optimal_k === nothing
            else
                @test result.optimal_k == r_profile.optimal_k[1]
                @test result.inc_threshold == r_profile.inc_threshold[1]
                @test result.ev_threshold == r_profile.ev_threshold[1]
            end
            for k in 1:n
                r_cut = cuts[(cuts.case .== case) .& (cuts.method .== r_method) .& (cuts.k .== k), :]
                @test cutree(context; k) == r_cut.membership
                @test cutree(context; k) == region_reference_cut(hcat(r_merge.left, r_merge.right), n, k)
            end
        end
    end
end

@testset "Prepared regionalization and diagnostics keep packed storage" begin
    n = 1000
    rng = MersenneTwister(201)
    distances = Matrix(Symmetric(rand(rng, n, n)))
    distances[diagind(distances)] .= 0
    context = prepare_phyloregion(distances)
    @test Base.summarysize(context) < 4_300_000
    cutree(context; k = 10)
    @test (@allocated cutree(context; k = 10)) < 100_000
    optimal_phyloregion(context)
    @test (@allocated optimal_phyloregion(context)) < 400_000
    select_linkage(context; methods = (:average,))
    @test (@allocated select_linkage(context; methods = (:average,))) < 400_000
end
