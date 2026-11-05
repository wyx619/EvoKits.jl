using Random, SparseArrays, Statistics, LinearAlgebra

# Independent small-tree oracle: original-edge paths and explicit tip pairs.
function species_reference(community, tree; scaled, precision = 256)
    return setprecision(BigFloat, precision) do
        edge_of = Dict(tree.child_of_edge[e] => e for e in 1:tree.nedges)
        paths = map(tree.tip_ids) do tip
            path = Set{Int}()
            while tip != tree.root
                edge = edge_of[tip]
                push!(path, edge)
                tip = tree.parent_of_node[tip]
            end
            path
        end
        heights = [sum((BigFloat(tree.edge_length[e]) for e in path); init = big"0") for path in paths]
        covariance, contrast = zeros(BigFloat, tree.ntips, tree.ntips), zeros(BigFloat, tree.ntips, tree.ntips)
        for i in 1:tree.ntips, j in 1:tree.ntips
            shared = sum((BigFloat(tree.edge_length[e]) for e in intersect(paths[i], paths[j])); init = big"0")
            covariance[i, j] = scaled ? shared / sqrt(heights[i] * heights[j]) : shared
            contrast[i, j] = scaled ? 1 - covariance[i, j] :
                sum((BigFloat(tree.edge_length[e]) for e in symdiff(paths[i], paths[j])); init = big"0") / 2
        end
        evenness, clustering, variability = [fill(NaN, size(community, 1)) for _ in 1:3]
        for row in axes(community, 1)
            tips = findall(>(0), community[row, :])
            m = length(tips)
            m > 1 || continue
            weights = BigFloat.(community[row, tips])
            weights ./= sum(weights)
            weighted = sum(weights[i] * weights[j] * contrast[tips[i], tips[j]]
                for j in 1:m for i in 1:(j - 1))
            evenness[row] = Float64(2weighted / (1 - big(1) / m))
            clustering[row] = Float64(1 - mean(maximum(covariance[i, j] for j in tips if j != i) for i in tips))
            variability[row] = Float64(mean(contrast[i, j] for j in tips for i in tips if i < j))
        end
        return (; pse = evenness, psc = clustering, psv = variability)
    end
end

species_match(a, b; rtol = 3e-12, atol = 1e-14) =
    size(a) == size(b) && all(isapprox.(a, b; rtol, atol, nans = true))

function species_random_tree(rng, p)
    parts = ["sp_$i:" * string(rand(rng, [0.125, 0.5, 1.0, 3.0])) for i in 1:p]
    while length(parts) > 1
        shuffle!(rng, parts)
        children = splice!(parts, 1:rand(rng, 2:min(5, length(parts))))
        push!(parts, "(" * join(children, ",") * "):" * string(rand(rng, [0.0, 0.25, 2.0])))
    end
    return only(parts) * ";"
end

@testset "PSE/PSC explicit covariance oracle and shared PSD results" begin
    rng = MersenneTwister(20261106)
    mktempdir() do dir
        path = joinpath(dir, "species.nwk")
        for _ in 1:28
            p, n = rand(rng, 3:14), rand(rng, 4:9)
            write_tree(path, species_random_tree(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:9, n, p)
            matrix[rand(rng, n, p) .< 0.55] .= 0
            matrix[1, :] .= 0
            matrix[2, :] .= 0
            matrix[2, 1] = 1
            columns, rows = randperm(rng, p), randperm(rng, n)
            permuted = matrix[rows, columns]
            saved = copy(permuted)
            edge_lengths = copy(tree.edge_length)
            for scaled in (false, true)
                reference = species_reference(matrix[rows, :], tree; scaled)
                for input in (permuted, sparse(permuted), view(permuted, :, :))
                    options = (; scale_vcv = scaled, tip_labels = tree.tip_labels[columns])
                    e, c = pse(input, tree; options...), psc(input, tree; options...)
                    v, r = psv(input, tree; options...), psr(input, tree; options...)
                    result = psd(input, tree; options...)
                    @test species_match(e.pse, reference.pse)
                    @test species_match(c.psc, reference.psc)
                    @test species_match(result.psv, reference.psv)
                    @test result.richness == e.richness == c.richness == vec(sum(input .> 0; dims = 2))
                    @test propertynames(result) == (:psv, :var_psv, :psc, :psr, :var_psr, :pse, :richness)
                    for (field, expected) in ((:psv, v.psv), (:var_psv, v.variance), (:psc, c.psc),
                        (:psr, r.psr), (:var_psr, r.variance), (:pse, e.pse))
                        @test species_match(getproperty(result, field), expected)
                        @test getproperty(result, field) isa Vector{Float64}
                    end
                    no_var = psd(input, tree; compute_var = false, options...)
                    @test all(isnan, no_var.var_psv) && all(isnan, no_var.var_psr)
                    for field in (:psv, :psc, :psr, :pse, :richness)
                        @test species_match(getproperty(result, field), getproperty(no_var, field))
                    end
                    data = CommunityData(input; species_labels = tree.tip_labels[columns])
                    labeled = psd(data, tree; scale_vcv = scaled)
                    @test all(species_match(getproperty(result, field), getproperty(labeled, field)) for field in propertynames(result))
                end
                presence = matrix .> 0
                @test species_match(pse(presence, tree; scale_vcv = scaled).pse,
                    psv(presence, tree; scale_vcv = scaled).psv)
            end
            @test permuted == saved
            @test tree.edge_length == edge_lengths
        end
    end
end

@testset "Species diversity extreme scales, zero paths and PSC definition" begin
    mktempdir() do dir
        path = joinpath(dir, "extreme.nwk")
        for (newick, matrix) in (
            ("(A:1,B:1);", [1e308 1e308; 1e-300 1e-300; 1e300 1e-300]),
            ("(A:1e300,B:1e300);", [1e300 1e-300; 1 nextfloat(0.0)]),
            ("((A:1e-100,B:1e-100):1e200,C:1e200);", [1 2 0; 2 1 0]),
            ("(((A:1,B:10000):1,C:0.01):1,D:1);", [1 1 1 1; 1 1 1 0]),
            ("((A:0,B:0):1,C:1);", [1 2 0; 1 1 1]),
            ("((A:1,B:1,C:1,D:1):1e308,E:1);", [1 1 1 1 0]))
            write_tree(path, newick)
            tree = read_tree(path)
            for scaled in (false, true)
                reference = species_reference(matrix, tree; scaled, precision = 2304)
                @test species_match(pse(sparse(matrix), tree; scale_vcv = scaled).pse, reference.pse; atol = 0, rtol = 2e-12)
                @test species_match(psc(matrix, tree; scale_vcv = scaled).psc, reference.psc; atol = 0, rtol = 2e-12)
                result = psd(matrix, tree; scale_vcv = scaled, compute_var = false)
                @test species_match(result.pse, reference.pse; atol = 0, rtol = 2e-12)
                @test species_match(result.psc, reference.psc; atol = 0, rtol = 2e-12)
            end
        end
        write_tree(path, "(A:1,B:1);")
        tree = read_tree(path)
        for matrix in (Float16[1000 2000], Float32[1 2], Real[1 2//1],
            [big(10)^1000 2big(10)^1000], [big"1e1000" big"2e1000"],
            Rational{BigInt}[big(10)^1000//1 2big(10)^1000//1])
            @test only(pse(matrix, tree).pse) ≈ 8 / 9
            @test only(psd(matrix, tree).pse) ≈ 8 / 9
            @test only(psc(matrix, tree).psc) == 1
        end
        # Values below Float64's weight range can collectively remain visible.
        native = reshape(BigFloat[1, big"1e-324", big"1e-324"], 1, :)
        write_tree(path, "(A:1,B:1,C:1);")
        tree = read_tree(path)
        expected = species_reference(native, tree; scaled = true).pse
        @test only(expected) > 0
        @test pse(native, tree).pse == expected
        result = psd([1 2 3], tree)
        arrays = [getproperty(result, field) for field in propertynames(result)]
        @test all(arrays[i] !== arrays[j] for j in eachindex(arrays) for i in 1:(j - 1))
    end
end

@testset "Species diversity validation, labels, empty inputs and sparse storage" begin
    mktempdir() do dir
        path = joinpath(dir, "species.nwk")
        write_tree(path, "((A:1,B:1):1,C:2);")
        tree = read_tree(path)
        explicit = sparse([1, 1, 2], [1, 2, 3], [0.0, 2.0, 1.0], 3, 3)
        saved = copy(explicit)
        for f in (pse, psc, psd)
            @test f(explicit, tree).richness == [1, 1, 0]
            @test all(isnan, getproperty(f(explicit, tree), f === psc ? :psc : :pse))
            @test isempty(getproperty(f(spzeros(0, 3), tree), f === psc ? :psc : :pse))
            data = CommunityData(sparse([1 1]); species_labels = ["C", "A"])
            @test f(data, tree).richness == [2]
            @test_throws ArgumentError f(data, tree; tip_labels = ["C", "A"])
            @test_throws ArgumentError f(ones(1, 2), tree)
            @test_throws ArgumentError f(ones(1, 3), tree; tip_labels = ["A", "A", "C"])
            @test_throws ArgumentError f(ones(1, 3), tree; tip_labels = ["A", "B", "X"])
            for bad in (-1.0, NaN, Inf)
                @test_throws ArgumentError f([bad 0 0], tree)
                invalid = deepcopy(tree)
                invalid.edge_length[1] = bad
                @test_throws ArgumentError f(ones(1, 3), invalid)
            end
            invalid = deepcopy(tree)
            invalid.edge_length .= 1e308
            @test_throws ArgumentError f(ones(1, 3), invalid)
        end
        @test explicit == saved && nonzeros(explicit) == nonzeros(saved)
        for newick in ("A;", "(A:0,B:0);")
            write_tree(path, newick)
            tree = read_tree(path)
            for f in (pse, psc, psd)
                @test_throws ArgumentError f(ones(1, tree.ntips), tree)
                @test f(ones(1, tree.ntips), tree; scale_vcv = false).richness == [tree.ntips]
            end
        end
        @test pse(ones(1, 2), tree; scale_vcv = false).pse == [0]
        @test psc(ones(1, 2), tree; scale_vcv = false).psc == [1]
        p, n = 20_000, 8
        write_tree(path, "(" * join(("sp_$i:1" for i in 1:p), ",") * ");")
        tree = read_tree(path)
        matrix = sparse(repeat(1:n, inner = 3), 1:3n, repeat([1, 2, 3], n), n, p)
        for scaled in (false, true)
            for (f, field) in ((pse, :pse), (psc, :psc), (psd, :pse))
                result = f(matrix, tree; scale_vcv = scaled)
                @test getproperty(result, field) ≈ fill(f === psc ? 1 : 11 / 12, n)
                @test @allocated(f(matrix, tree; scale_vcv = scaled)) < 12_000_000
            end
            combined = @allocated psd(matrix, tree; scale_vcv = scaled)
            separate = sum(@allocated(f(matrix, tree; scale_vcv = scaled)) for f in (psv, psr, pse, psc))
            @test combined < separate
        end
    end
end
