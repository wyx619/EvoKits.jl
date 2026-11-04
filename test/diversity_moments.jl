using SparseArrays, Random, LinearAlgebra, Statistics

# Independent original-edge root paths, used only on small trees.
function moment_paths(tree)
    edge_of = Dict(tree.child_of_edge[e] => e for e in 1:tree.nedges)
    return map(tree.tip_ids) do tip
        path = Set{Int}()
        while tip != tree.root
            edge = edge_of[tip]
            push!(path, edge)
            tip = tree.parent_of_node[tip]
        end
        path
    end
end

function moment_newick(rng, p)
    parts = ["sp_$i:" * string(rand(rng, [0.125, 0.5, 1.0, 3.0])) for i in 1:p]
    while length(parts) > 1
        shuffle!(rng, parts)
        children = splice!(parts, 1:rand(rng, 2:min(4, length(parts))))
        push!(parts, "(" * join(children, ",") * "):" * string(rand(rng, [0.0, 0.25, 2.0])))
    end
    return only(parts) * ";"
end

function psv_enumeration(tree; scaled = true)
    paths = moment_paths(tree)
    heights = [sum((BigFloat(tree.edge_length[e]) for e in path); init = BigFloat(0)) for path in paths]
    p = tree.ntips
    contrast = zeros(BigFloat, p, p)
    for j in 1:p, i in 1:(j - 1)
        shared = sum((BigFloat(tree.edge_length[e]) for e in intersect(paths[i], paths[j])); init = BigFloat(0))
        contrast[i, j] = contrast[j, i] = scaled ?
            1 - shared / sqrt(heights[i] * heights[j]) :
            sum((BigFloat(tree.edge_length[e]) for e in symdiff(paths[i], paths[j])); init = BigFloat(0)) / 2
    end
    rows = zeros(Int, 2^p, p)
    values = fill(BigFloat(NaN), 2^p)
    by_richness = [BigFloat[] for _ in 1:p]
    for mask in 0:(2^p - 1)
        tips = findall(i -> !iszero(mask & (1 << (i - 1))), 1:p)
        rows[mask + 1, tips] .= 1
        length(tips) > 1 || continue
        pairs = [contrast[i, j] for j in tips for i in tips if i < j]
        values[mask + 1] = mean(pairs)
        push!(by_richness[length(tips)], values[mask + 1])
    end
    variance = [isempty(v) ? NaN : Float64(mean((v .- mean(v)) .^ 2)) for v in by_richness]
    return rows, Float64.(values), variance
end

# Enumerate all independent Bernoulli tip outcomes; no ancestor recurrence.
function expected_pd_enumeration(tree)
    paths = moment_paths(tree)
    p = tree.ntips
    outcomes = map(0:(2^p - 1)) do mask
        tips = findall(i -> !iszero(mask & (1 << (i - 1))), 1:p)
        edges = union(Set{Int}(), (paths[i] for i in tips)...)
        value = sum((BigFloat(tree.edge_length[e]) for e in edges); init = BigFloat(0))
        (; k = length(tips), value)
    end
    means, variances = BigFloat[], BigFloat[]
    for m in 1:p
        probability = BigFloat(m) / p
        weights = [probability^o.k * (1 - probability)^(p - o.k) for o in outcomes]
        center = sum(w * o.value for (w, o) in zip(weights, outcomes))
        push!(means, center)
        push!(variances, sum(w * (o.value - center)^2 for (w, o) in zip(weights, outcomes)))
    end
    return Float64.(means), Float64.(variances)
end

@testset "PSV/PSR exact finite-pool variance and leaf-order invariance" begin
    rng = MersenneTwister(20261105)
    mktempdir() do dir
        path = joinpath(dir, "moments.nwk")
        for trial in 1:22
            p = rand(rng, 2:8)
            write_tree(path, moment_newick(rng, p))
            tree = read_tree(path)
            for scaled in (false, true)
                rows, expected, variance = psv_enumeration(tree; scaled)
                richness = vec(sum(rows; dims = 2))
                expected_var = [k > 1 ? variance[k] : NaN for k in richness]
                columns = randperm(rng, p)
                for input in (rows, sparse(rows[:, columns]))
                    labels = issparse(input) ? tree.tip_labels[columns] : tree.tip_labels
                    actual = psv(input, tree; scale_vcv = scaled, tip_labels = labels)
                    @test all(isapprox.(actual.psv, expected; atol = 1e-14, rtol = 2e-12, nans = true))
                    @test all(isapprox.(actual.variance, expected_var; atol = 1e-14, rtol = 2e-12, nans = true))
                    @test actual.richness == richness
                    @test actual.variance[end] == 0
                    @test all(x -> isnan(x) || x >= 0, actual.variance)
                    r = psr(input, tree; scale_vcv = scaled, tip_labels = labels)
                    @test isequal(r.psr, actual.psv .* richness)
                    @test isequal(r.variance, actual.variance .* richness .^ 2)
                    @test all(isnan, psv(input, tree; scale_vcv = scaled,
                        compute_var = false, tip_labels = labels).variance)
                end
            end
        end
        results = []
        for newick in ("((A:1,B:1):1,C:2);", "(C:2,(B:1,A:1):1);")
            write_tree(path, newick)
            tree = read_tree(path)
            push!(results, psv([1 1 0; 1 1 1], tree))
            @test results[end].variance[1] ≈ 1 / 18
            @test results[end].variance[2] == 0
        end
        @test results[1].variance ≈ results[2].variance
    end
end

@testset "PSV sparse input, empty sites, validation and short branches" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "((A:1,B:1):1,C:2);")
        tree = read_tree(path)
        matrix = [1 1 0; 0 0 0; 0 0 3]
        saved = copy(matrix)
        actual = psv(view(matrix, :, :), tree)
        @test actual.psv[1] ≈ 0.5
        @test all(isnan, actual.psv[2:3])
        @test matrix == saved
        data = CommunityData(sparse(matrix[:, [3, 1]]); species_labels = ["C", "A"])
        @test isequal(psv(data, tree).richness, [1, 0, 1])
        @test isempty(psv(zeros(0, 3), tree).psv)
        explicit = sparse([1, 1], [1, 2], [0.0, 1.0], 1, 3)
        @test psv(explicit, tree).richness == [1]
        @test nnz(explicit) == 2
        for bad in (-1.0, Inf, NaN)
            @test_throws ArgumentError psv([bad 0 0], tree)
            invalid = deepcopy(tree)
            invalid.edge_length[1] = bad
            @test_throws ArgumentError psv(ones(1, 3), invalid)
        end
        @test_throws ArgumentError psv(ones(1, 2), tree)
        @test_throws ArgumentError psv(ones(1, 3), tree; tip_labels = ["A", "A", "C"])
        write_tree(path, "(A:0,B:0);")
        tree = read_tree(path)
        @test_throws ArgumentError psv(ones(1, 2), tree)
        @test psv(ones(1, 2), tree; scale_vcv = false).psv == [0]
        write_tree(path, "A;")
        tree = read_tree(path)
        @test isnan(only(psv(ones(1, 1), tree; scale_vcv = false).variance))
        write_tree(path, "((A:1e-100,B:1e-100):1e200,C:1e200);")
        tree = read_tree(path)
        for scaled in (false, true)
            @test psv(sparse([1 1 0]), tree; scale_vcv = scaled).psv[1] ≈ (scaled ? 1e-300 : 1e-100) rtol=1e-12
        end
        p, n = 20_000, 8
        write_tree(path, "(" * join(("sp_$i:1" for i in 1:p), ",") * ");")
        tree = read_tree(path)
        matrix = sparse(repeat(1:n, inner = 3), 1:3n, ones(3n), n, p)
        for scaled in (false, true)
            result = psv(matrix, tree; scale_vcv = scaled)
            @test result.psv == ones(n)
            @test result.variance == zeros(n)
            @test @allocated(psv(matrix, tree; scale_vcv = scaled)) < 12_000_000
        end
    end
end

@testset "Expected PD exact Bernoulli outcomes and upper bound" begin
    rng = MersenneTwister(20261106)
    mktempdir() do dir
        path = joinpath(dir, "expected.nwk")
        for _ in 1:24
            p = rand(rng, 2:8)
            write_tree(path, moment_newick(rng, p))
            tree = read_tree(path)
            saved = copy(tree.edge_length)
            means, variances = expected_pd_enumeration(tree)
            exact = expected_pd(tree; upper_bound = false)
            bound = expected_pd(tree)
            @test exact.expected_pd ≈ means rtol=2e-13
            @test exact.variance_pd ≈ variances rtol=2e-13
            @test bound.expected_pd ≈ means rtol=2e-13
            @test all(bound.variance_pd .>= exact.variance_pd .- 1e-12)
            @test exact.variance_pd[end] == bound.variance_pd[end] == 0
            @test exact.richness == collect(1:p)
            @test tree.edge_length == saved
            @test pd_variance(tree; upper_bound = false).variance_pd == exact.variance_pd
            @test ead(tree).num_children == edge_abundance_distribution(tree).num_children
        end
        for newick in ("A;", "(A:0,B:0,C:0);", "(A:1,B:1,C:1);")
            write_tree(path, newick)
            tree = read_tree(path)
            means, variances = expected_pd_enumeration(tree)
            for upper_bound in (false, true)
                actual = expected_pd(tree; upper_bound)
                @test actual.expected_pd ≈ means
                @test actual.variance_pd ≈ variances
            end
        end
        write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
        tree = read_tree(path)
        @test expected_pd(tree).expected_pd ≈ [1.875, 3.5, 4.875, 6]
        @test expected_pd(tree).variance_pd ≈ [3.2109375, 2.875, 1.3359375, 0]
        @test expected_pd(tree; upper_bound = false).variance_pd ≈ [2.3671875, 2.375, 1.2421875, 0]
        for bad in (-1.0, Inf, NaN)
            invalid = deepcopy(tree)
            invalid.edge_length[1] = bad
            @test_throws ArgumentError expected_pd(invalid)
        end
        p = 4_000
        left = "(" * join(("sp_$i:1e-100" for i in 1:(p ÷ 2)), ",") * "):1e300"
        right = join(("sp_$i:1e-100" for i in (p ÷ 2 + 1):p), ",")
        write_tree(path, "(" * left * "," * right * ");")
        tree = read_tree(path)
        for upper_bound in (false, true)
            result = expected_pd(tree; upper_bound)
            @test !any(isnan, result.variance_pd)
            @test isfinite(result.variance_pd[p ÷ 2])
            @test result.variance_pd[p ÷ 2] ≈ Float64(big"1e600" * big"0.5"^2000) rtol=1e-12
            @test result.variance_pd[end] == 0
            @test @allocated(expected_pd(tree; upper_bound)) < 3_000_000
        end
    end
end
