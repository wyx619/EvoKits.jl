# Reuses the independent original-edge and subset-enumeration oracles from
# phylogenetic_species.jl/diversity_moments.jl, plus the iterative comb fixture.
function psv_virtual_full(matrix, tree, scaled, compute_var)
    rows = EvoKits._community_tip_rows(matrix, tree, nothing)
    geometry = EvoKits._psv_geometry(tree, scaled)
    return EvoKits._psv_values(rows, tree, geometry, zeros(Int, tree.nnodes),
        zeros(tree.nnodes), zeros(tree.nnodes); compute_var)
end

@testset "PSV connecting trees retain original richness and exact full-pool variance" begin
    mktempdir() do dir
        path = joinpath(dir, "psv.nwk")
        for newick in ("(((A:1,B:2):0,C:4):3,(D:3,E:7,F:11):2);",
            "((((A:1):2,B:3):4,C:2):1,D:2,E:5,F:8);")
            write(path, newick)
            tree = read_tree(path)
            # E is observed only in a singleton; F is never observed. Both
            # still belong to the finite population for variance.
            matrix = zeros(Int, 18, 6)
            matrix[2, 5] = 9
            for mask in 0:15, tip in 1:4
                matrix[mask + 3, tip] = isodd(mask >> (tip - 1)) ? tip : 0
            end
            richness = vec(sum(matrix .> 0; dims = 2))
            for scaled in (false, true)
                _, all_values, pool_variance = psv_enumeration(tree; scaled)
                expected = [all_values[1 + sum((1 << (tip - 1) for tip in findall(>(0), row)); init = 0)]
                    for row in eachrow(matrix)]
                expected_var = [r > 1 ? pool_variance[r] : NaN for r in richness]
                rows = EvoKits._community_tip_rows(sparse(matrix), tree, nothing)
                geometry = EvoKits._psv_geometry(tree, scaled)
                for strategy in (:full, :virtual, :auto), compute_var in (false, true)
                    plan = EvoKits._species_plan(rows, tree, geometry; strategy)
                    actual = EvoKits._species_psv_values(rows, plan, compute_var ? pool_variance : Float64[])
                    @test species_match(actual.psv, expected)
                    @test species_match(actual.variance, compute_var ? expected_var : fill(NaN, 18))
                    @test actual.richness == richness
                    @test actual.scale_vcv == scaled
                    strategy === :virtual && @test isempty(plan.rows[2])
                end
                order = [6, 2, 5, 1, 4, 3]
                for input in (matrix[:, order], sparse(matrix[:, order]), view(matrix, :, order))
                    actual = psv(input, tree; scale_vcv = scaled, tip_labels = tree.tip_labels[order])
                    @test species_match(actual.psv, expected)
                    @test species_match(actual.variance, expected_var)
                    @test actual.richness == richness
                end
                data = CommunityData(sparse(matrix[:, 1:5]); species_labels = tree.tip_labels[1:5])
                @test species_match(psv(data, tree; scale_vcv = scaled).variance, expected_var)
            end
        end
    end
end

@testset "PSV automatic sparse and mixed site paths match original-tree moments" begin
    tree = species_comb(4096)
    n = 64
    matrix = sparse(repeat(1:n; inner = 3), [1 + mod(719i, tree.ntips) for i in 1:(3n)],
        repeat([1, 2, 3], n), n, tree.ntips)
    matrix[1, :] .= 0
    matrix[2, :] .= 0
    matrix[2, end] = 1
    saved = copy(matrix)
    for mixed in (false, true)
        input = mixed ? vcat(ones(Int, 1, tree.ntips), matrix, ones(Int, 1, tree.ntips), matrix) : matrix
        for scaled in (false, true), compute_var in (false, true)
            rows = EvoKits._community_tip_rows(input, tree, nothing)
            plan = EvoKits._species_plan(rows, tree, EvoKits._psv_geometry(tree, scaled))
            @test plan.index !== nothing
            @test mixed ? plan.tree.nnodes == tree.nnodes : plan.tree.nnodes < tree.nnodes ÷ 4
            expected = psv_virtual_full(input, tree, scaled, compute_var)
            actual = psv(input, tree; scale_vcv = scaled, compute_var)
            @test species_match(actual.psv, expected.psv)
            @test isequal(actual.variance, expected.variance)
            @test actual.richness == expected.richness
            r = psr(input, tree; scale_vcv = scaled, compute_var)
            @test species_match(r.psr, expected.psv .* expected.richness)
            @test isequal(r.variance, expected.variance .* expected.richness .^ 2)
            combined = psd(input, tree; scale_vcv = scaled, compute_var)
            @test species_match(actual.psv, combined.psv)
            @test isequal(actual.variance, combined.var_psv)
            @test actual.psv !== actual.variance
        end
    end
    @test matrix == saved
    # Every species occurs somewhere: the global tree cannot be compressed,
    # but each two-tip query still has just three connecting nodes.
    tree = species_comb(1024)
    matrix = sparse(repeat(1:512; inner = 2), 1:1024, ones(1024), 512, 1024)
    for scaled in (false, true)
        actual = psv(matrix, tree; scale_vcv = scaled)
        for site in 1:512
            a, b = 2site - 1, 2site
            expected = scaled ? 1 - (a - 1) / sqrt(a * (b == tree.ntips ? b - 1 : b)) :
                (b == tree.ntips ? 2 : b - a + 2) / 2
            @test actual.psv[site] ≈ expected rtol = 5e-12
        end
        @test length(unique(actual.variance)) == 1
    end
end

@testset "PSV connecting paths preserve extreme shared stems and input validation" begin
    mktempdir() do dir
        path = joinpath(dir, "psv.nwk")
        for newick in ("((A:1e-100,B:2e-100):1e200,C:1e200);",
            "((((A:0,B:0):0):1,C:2):3,D:4);")
            write(path, newick)
            tree = read_tree(path)
            matrix = zeros(Int, 3, tree.ntips)
            matrix[1, 1:2] .= 1
            matrix[3, end] = 1
            for scaled in (false, true)
                reference = species_reference(matrix, tree; scaled, precision = 2304).psv
                rows = EvoKits._community_tip_rows(matrix, tree, nothing)
                plan = EvoKits._species_plan(rows, tree, EvoKits._psv_geometry(tree, scaled); strategy = :virtual)
                actual = EvoKits._species_psv_values(rows, plan, Float64[])
                @test species_match(actual.psv, reference; atol = 0, rtol = 3e-12)
                @test actual.richness == [2, 0, 1]
            end
        end
        write(path, "((A:1,B:2):1,C:0);")
        tree = read_tree(path)
        @test_throws ArgumentError psv([1 1 0], tree)
        @test all(isfinite, psv([1 1 0], tree; scale_vcv = false).psv)
    end
end

@testset "PSV connecting storage on a deep sparse tree" begin
    tree = species_comb(20_000)
    matrix = sparse(repeat(1:64; inner = 3), [1 + mod(719i, tree.ntips) for i in 1:192],
        ones(192), 64, tree.ntips)
    for scaled in (false, true), compute_var in (false, true)
        psv(matrix, tree; scale_vcv = scaled, compute_var)
        @test @allocated(psv(matrix, tree; scale_vcv = scaled, compute_var)) < 12_000_000
    end
    rows = EvoKits._community_tip_rows(matrix, tree, nothing)
    plan = EvoKits._species_plan(rows, tree, EvoKits._psv_geometry(tree, true))
    @test plan.tree.nnodes <= 2nnz(matrix) - 1
    @test Base.summarysize(plan) < 300_000
    workspace = EvoKits._SpeciesSiteWorkspace(plan.tree.nnodes)
    @test all(plan.rows) do tips
        local_tree, _, _ = EvoKits._species_site!(workspace, plan, tips)
        local_tree.nnodes <= 2length(tips) - 1
    end
end
