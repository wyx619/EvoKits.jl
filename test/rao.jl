using SparseArrays, Random, LinearAlgebra

# Independent definition: enumerate tip pairs using original edges on the
# symmetric difference of their root paths, with BigFloat abundance arithmetic.
# No production distance kernel, root-depth subtraction or LCA is used.
function rao_oracle(matrix, tree; weighted = true)
    paths = [Set{Int}() for _ in tree.tip_ids]
    edge_of = Dict(tree.child_of_edge[e] => e for e in 1:tree.nedges)
    for (node, path) in zip(tree.tip_ids, paths)
        while node != tree.root
            edge = edge_of[node]
            push!(path, edge)
            node = tree.parent_of_node[node]
        end
    end
    distances = [sum((BigFloat(tree.edge_length[e]) for e in symdiff(a, b)); init = BigFloat(0)) / 2
        for a in paths, b in paths]
    raw = BigFloat.(Matrix(matrix))
    site_totals = vec(sum(raw; dims = 2))
    composition = weighted ? copy(raw) : BigFloat.(raw .> 0)
    composition ./= sum(composition; dims = 2)
    n, p = size(raw)
    mean_distance(a, b) = sum((a[i] * b[j] * distances[i, j] for i in 1:p, j in 1:p);
        init = BigFloat(0))
    dkl = [mean_distance(composition[i, :], composition[j, :]) for i in 1:n, j in 1:n]
    dkk = diag(dkl)
    h = [dkl[i, j] - (dkk[i] + dkk[j]) / 2 for i in 1:n, j in 1:n]
    pooled = vec(sum(raw; dims = 1)) ./ sum(raw)
    total = mean_distance(pooled, pooled)
    alpha = sum(dkk .* (site_totals ./ sum(raw)))
    beta = total - alpha
    return RaoDResult(Dkk = Float64.(dkk), Dkl = Float64.(dkl), H = Float64.(h),
        total = Float64(total), alpha = Float64(alpha), beta = Float64(beta), Fst = Float64(beta / total))
end

function test_rao_result(actual, expected)
    tol = 3e-12
    distance_scale = maximum((abs(x) for x in expected.Dkl if isfinite(x)); init = 0.0)
    partition_scale = max(isfinite(expected.total) ? abs(expected.total) : 0,
        isfinite(expected.alpha) ? abs(expected.alpha) : 0)
    for field in (:Dkk, :Dkl, :H, :total, :alpha, :beta, :Fst)
        a, b = getproperty(actual, field), getproperty(expected, field)
        atol = field === :H ? tol * distance_scale :
            field === :beta ? tol * partition_scale :
            field === :Fst && expected.total > 0 ? tol * partition_scale / expected.total : 0.0
        @test isapprox(a, b; atol, rtol = tol, nans = true)
    end
end

function rao_test_newick(rng, p)
    parts = [("sp_$i", 0.0) for i in 1:p]
    while length(parts) > 1
        shuffle!(rng, parts)
        children = splice!(parts, 1:rand(rng, 2:min(4, length(parts))))
        height = maximum(last, children) + rand(rng, [0.0, 0.125, 0.5, 1.0, 3.0])
        parent = "(" * join((s * ":" * string(height - h) for (s, h) in children), ",") * ")"
        push!(parts, (parent, height))
    end
    return first(only(parts)) * ";"
end

@testset "RaoD analytic partition and original-abundance pooling" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "(A:1,B:1,C:1);")
        tree = read_tree(path)
        matrix = [1 1 0; 0 1 1]
        for weighted in (false, true), input in (matrix, sparse(matrix), view(matrix, :, :))
            actual = raoD(input, tree; abundance_weighted = weighted)
            test_rao_result(actual, RaoDResult(Dkk = [0.5, 0.5], Dkl = [0.5 0.75; 0.75 0.5],
                H = [0.0 0.25; 0.25 0.0], total = 0.625, alpha = 0.5, beta = 0.125, Fst = 0.2))
            @test actual.Dkk == diag(actual.Dkl)
            @test isequal(actual.Dkl, transpose(actual.Dkl))
            @test isequal(actual.H, transpose(actual.H))
        end
        skewed = [1000 1 0; 1000 1 0]
        incidence = raoD(sparse(skewed), tree; abundance_weighted = false)
        abundance = raoD(skewed, tree)
        @test incidence.alpha == 0.5
        @test incidence.total ≈ 2000 / 1001^2
        @test incidence.total ≈ abundance.total
        @test incidence.beta < 0 && incidence.Fst < 0
        test_rao_result(incidence, rao_oracle(skewed, tree; weighted = false))
        # Returned matrices/vectors own their storage independently.
        actual = raoD(matrix, tree)
        actual.Dkl[1, 1] = -1
        @test actual.Dkk[1] == 0.5 && actual.H[1, 1] == 0
        actual.H[1, 2] = -2
        @test actual.Dkl[1, 2] == 0.75
    end
end

@testset "RaoD empty, singleton and zero-diversity conventions" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        for newick in ("(A:1,B:1,C:1);", "(A:0,B:0,C:0);", "A;")
            write_tree(path, newick)
            tree = read_tree(path)
            p = tree.ntips
            mixed = zeros(3, p)
            mixed[1, :] .= 1
            mixed[2, 1] = 2
            for weighted in (false, true), input in (mixed, sparse(mixed), spzeros(2, p), zeros(0, p))
                actual = raoD(input, tree; abundance_weighted = weighted)
                test_rao_result(actual, rao_oracle(input, tree; weighted))
                @test size(actual.Dkl) == (size(input, 1), size(input, 1))
                if size(input, 1) == 3
                    @test actual.Dkk[2] == 0
                    @test all(isnan, actual.Dkl[3, :]) && all(isnan, actual.H[:, 3])
                    @test isfinite(actual.total)
                    @test all(isnan, (actual.alpha, actual.beta, actual.Fst))
                elseif size(input, 1) == 0
                    @test isempty(actual.Dkk) && isempty(actual.H)
                    @test actual.alpha == 0 && isnan(actual.total)
                else
                    @test all(isnan, actual.Dkk)
                    @test all(isnan, (actual.total, actual.alpha, actual.beta, actual.Fst))
                end
            end
            mono = zeros(2, p)
            mono[:, 1] .= [1, 3]
            actual = raoD(sparse(mono), tree)
            @test actual.total == actual.alpha == actual.beta == 0
            @test isnan(actual.Fst)
        end
        write_tree(path, "(A:1,B:1,C:1);")
        tree = read_tree(path)
        explicit = sparse([1, 1, 2], [1, 2, 3], [0.0, 2, 3], 3, 3)
        test_rao_result(raoD(explicit, tree), rao_oracle(explicit, tree))
        @test nnz(explicit) == 3
    end
end

@testset "RaoD randomized ultrametric trees, input storage and alignment" begin
    rng = MersenneTwister(20261108)
    mktempdir() do dir
        path = joinpath(dir, "random.nwk")
        for trial in 1:24
            n, p = rand(rng, 2:9), rand(rng, 3:18)
            write_tree(path, rao_test_newick(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:8, n, p) ./ 4
            matrix[rand(rng, n, p) .< 0.65] .= 0
            iseven(trial) && (matrix[end, :] .= 0)
            columns, sites = randperm(rng, p), randperm(rng, n)
            saved_lengths = copy(tree.edge_length)
            for weighted in (false, true)
                expected = rao_oracle(matrix, tree; weighted)
                for (input, tips) in ((matrix, tree.tip_labels), (view(matrix, :, :), tree.tip_labels),
                    (sparse(matrix), tree.tip_labels), (sparse(matrix[:, columns]), tree.tip_labels[columns]))
                    saved = copy(input)
                    test_rao_result(raoD(input, tree; tip_labels = tips, abundance_weighted = weighted), expected)
                    @test input == saved
                end
                data = CommunityData(sparse(matrix[sites, columns]); species_labels = tree.tip_labels[columns])
                test_rao_result(raoD(data, tree; abundance_weighted = weighted),
                    rao_oracle(matrix[sites, :], tree; weighted))
            end
            @test tree.edge_length == saved_lengths
        end
        write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
        tree = read_tree(path)
        data = CommunityData(sparse([2 1; 0 3]); species_labels = ["C", "A"])
        test_rao_result(raoD(data, tree), rao_oracle([1 0 2 0; 3 0 0 0], tree))
    end
end

@testset "RaoD numeric types, abundance scaling and short descendant paths" begin
    mktempdir() do dir
        path = joinpath(dir, "types.nwk")
        write_tree(path, "((A:1,B:1):1,C:2);")
        tree = read_tree(path)
        for matrix in (Float16[1 3 0; 0 2 7], Float32[0.1 0.7 0; 0.3 0 0.9],
            BigInt[big(10)^400 big(3)*big(10)^400 0; 0 2 7],
            BigFloat[big"1e400" big"3e400" 0; 0 big"1e-400" big"7e-400"],
            Rational{Int}[1//3 2//3 0; 0 1//7 6//7], Bool[1 1 0; 0 1 1],
            Real[1 3 0; 0 2 7])
            for input in (matrix, sparse(matrix)), weighted in (false, true)
                test_rao_result(raoD(input, tree; abundance_weighted = weighted),
                    rao_oracle(input, tree; weighted))
            end
        end
        matrix = [1.0 3 0; 0 2 7]
        for scale in (1e-300, 1e300), weighted in (false, true)
            test_rao_result(raoD(sparse(matrix .* scale), tree; abundance_weighted = weighted),
                rao_oracle(matrix, tree; weighted))
        end
        # Large finite row totals need not themselves fit Float64.
        matrix = [1e308 1e308 0; 0 1e308 1e308]
        for weighted in (false, true)
            test_rao_result(raoD(sparse(matrix), tree; abundance_weighted = weighted),
                rao_oracle(matrix, tree; weighted))
        end
        # Cached root depths lose the terminal paths, but original branches do not.
        write_tree(path, "((A:1e-100,B:1e-100):1e200,C:1e200);")
        tree = read_tree(path)
        matrix = [2.0 1 0; 0 4 0; 1 0 0]
        for weighted in (false, true)
            actual = raoD(sparse(matrix), tree; abundance_weighted = weighted)
            test_rao_result(actual, rao_oracle(matrix, tree; weighted))
            @test actual.Dkk[1] > 0 && actual.total > 0
        end
        # Halve the final weighted sum, not subnormal weights before multiplying.
        write_tree(path, "(A:1e307,B:1e307);")
        tree = read_tree(path)
        matrix = [1.0 nextfloat(0.0)]
        actual = raoD(matrix, tree)
        test_rao_result(actual, rao_oracle(matrix, tree))
        @test actual.total > 0 && actual.beta == 0
    end
end

@testset "RaoD validation, ultrametric tolerance and sparse storage" begin
    mktempdir() do dir
        path = joinpath(dir, "validation.nwk")
        write_tree(path, "(A:1,B:1,C:1);")
        tree = read_tree(path)
        for bad in (-1.0, Inf, NaN), input in ([bad 0 0], sparse([bad 0 0]))
            @test_throws ArgumentError raoD(input, tree)
        end
        @test_throws ArgumentError raoD(ones(1, 2), tree)
        for labels in (["A", "B"], ["A", "B", "X"], ["A", "A", "C"])
            @test_throws ArgumentError raoD(ones(1, 3), tree; tip_labels = labels)
        end
        for bad in (-1.0, Inf, NaN)
            invalid = deepcopy(tree)
            invalid.edge_length[1] = bad
            @test_throws ArgumentError raoD(ones(1, 3), invalid)
        end
        invalid = deepcopy(tree)
        invalid.edge_length .= floatmax(Float64) / 2
        @test_throws ArgumentError raoD(ones(1, 3), invalid)
        write_tree(path, "(A:1,B:1.000000005);")
        @test isfinite(raoD(ones(1, 2), read_tree(path)).total)
        write_tree(path, "(A:1,B:1.00000002);")
        @test_throws ArgumentError raoD(ones(1, 2), read_tree(path))
        p, n = 20_000, 12
        write_tree(path, "(" * join(("sp_$i:1" for i in 1:p), ",") * ");")
        tree = read_tree(path)
        matrix = sparse(repeat(1:n, inner = 3), collect(1:3n), fill(1.0, 3n), n, p)
        for weighted in (false, true)
            actual = raoD(matrix, tree; abundance_weighted = weighted)
            @test actual.Dkk ≈ fill(2 / 3, n)
            @test actual.total ≈ 1 - 1 / (3n)
            @test actual.alpha ≈ 2 / 3
            @test @allocated(raoD(matrix, tree; abundance_weighted = weighted)) < 15_000_000
        end
    end
end
