# Uses the independent original-edge BigFloat oracle in phylogenetic_species.jl.
function species_prepared(matrix, tree, scaled; strategy = :auto)
    rows = EvoKits._species_abundance_rows(matrix, tree, nothing)
    geometry = EvoKits._psv_geometry(tree, scaled; need_heights = !scaled)
    plan = EvoKits._species_plan(rows.tips, tree, geometry; strategy)
    return rows, plan
end

# Iterative comb: exercise deep paths without recursion in the Newick parser.
function species_comb(p)
    n = 2p - 1
    parents, children = zeros(Int32, n), [Int32[] for _ in 1:n]
    for i in 1:(p - 1)
        node = Int32(p + i)
        children[node] = Int32[i, i == p - 1 ? p : node + 1]
        for child in children[node]
            parents[child] = node
        end
    end
    root = Int32(p + 1)
    preorder, depths, stack = Int32[], zeros(n), Int32[root]
    while !isempty(stack)
        node = pop!(stack)
        push!(preorder, node)
        for child in Iterators.reverse(children[node])
            depths[child] = depths[node] + 1
            push!(stack, child)
        end
    end
    pe, ce, first, last = Int32[], Int32[], zeros(Int32, n), zeros(Int32, n)
    for node in 1:n, child in children[node]
        push!(pe, node)
        push!(ce, child)
        first[node] == 0 && (first[node] = length(pe))
        last[node] = length(pe)
    end
    labels = ["sp_$i" for i in 1:p]
    tips = vcat(trues(p), falses(p - 1))
    return CompactTree(p, n, n - 1, pe, ce, ones(n - 1), root, parents,
        depths, tips, Int32.(1:p), reverse(preorder), Int32.(n:-1:(p + 1)),
        preorder, children, first, last, Dict(zip(labels, Int32.(1:p))),
        labels, vcat(labels, fill("", p - 1)))
end

@testset "Species connecting trees against independent covariance oracle" begin
    rng = MersenneTwister(20261107)
    mktempdir() do dir
        path = joinpath(dir, "virtual.nwk")
        for _ in 1:24
            p, n = rand(rng, 4:24), 7
            write_tree(path, species_random_tree(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:9, n, p)
            matrix[rand(rng, n, p) .< 0.8] .= 0
            matrix[1, :] .= 0
            matrix[2, :] .= 0
            matrix[2, p] = 1
            matrix[3, :] .= 1
            for scaled in (false, true)
                reference = species_reference(matrix, tree; scaled)
                variance_by_richness = psv([Int(j <= i) for i in 1:p, j in 1:p], tree;
                    scale_vcv = scaled).variance
                expected_variance = psv(matrix, tree; scale_vcv = scaled).variance
                for strategy in (:full, :virtual, :auto)
                    rows, plan = species_prepared(sparse(matrix), tree, scaled; strategy)
                    @test species_match(EvoKits._species_pse_values(rows, plan), reference.pse)
                    @test species_match(EvoKits._species_psc_values(plan), reference.psc)
                    combined = EvoKits._species_psd_values(rows, plan, variance_by_richness)
                    @test species_match(combined.psv, reference.psv)
                    @test species_match(combined.pse, reference.pse)
                    @test species_match(combined.psc, reference.psc)
                    @test species_match(combined.var_psv, expected_variance)
                    @test species_match(combined.psr, reference.psv .* combined.richness)
                    @test species_match(combined.var_psr, expected_variance .* combined.richness .^ 2)
                    # Changing leaf order changes traversal, not the statistic.
                    columns = randperm(rng, p)
                    options = (; scale_vcv = scaled, tip_labels = tree.tip_labels[columns])
                    @test species_match(pse(view(matrix, :, columns), tree; options...).pse, reference.pse)
                    @test species_match(psc(sparse(matrix[:, columns]), tree; options...).psc, reference.psc)
                end
            end
        end
        tree = species_comb(8)
        rows, plan = species_prepared(ones(1, 8), tree, true)
        @test_throws ArgumentError EvoKits._species_plan(rows.tips, tree, plan.geometry; strategy = :unknown)
    end
end

@testset "Species contraction preserves root stems, tiny contrasts and identity" begin
    mktempdir() do dir
        path = joinpath(dir, "stems.nwk")
        cases = (
            ("((((A:1e-100):2e-100,B:4e-100):1e200,C:1e200):1e200,D:1e200);", [1 2 0 0; 1 1 3 0]),
            ("(((A:1,B:3):0.125,C:7):1e16,D:1);", [1 3 0 0; 1 2 3 0]),
            ("((((A:1,B:2):1):1):4,C:2);", [1 3 0; 1 2 4]),
            ("(((A:0,B:0):3,C:1):4,D:2);", [1 2 0 0; 1 2 4 0]),
            ("((A:1e300,B:1e300):1e300,C:1);", [1e300 1e-300 0; 1 nextfloat(0.0) 0]),
            ("((A:1,B:1,C:1):1e308,D:1);", [1 1 1 0; 1 2 3 0]),
        )
        for (newick, matrix) in cases
            write_tree(path, newick)
            tree = read_tree(path)
            for scaled in (false, true), strategy in (:full, :virtual, :auto)
                reference = species_reference(matrix, tree; scaled, precision = 2304)
                rows, plan = species_prepared(sparse(matrix), tree, scaled; strategy)
                actual = EvoKits._species_pse_values(rows, plan)
                @test species_match(actual, reference.pse; atol = 0, rtol = 4e-12)
                @test species_match(EvoKits._species_psc_values(plan), reference.psc; atol = 0, rtol = 4e-12)
                combined = EvoKits._species_psd_values(rows, plan, Float64[])
                @test species_match(combined.pse, reference.pse; atol = 0, rtol = 4e-12)
                @test species_match(combined.psc, reference.psc; atol = 0, rtol = 4e-12)
                @test species_match(combined.psv, reference.psv; atol = 0, rtol = 4e-12)
            end
        end
        write_tree(path, "((A:1,B:2):4,C:0);")
        tree = read_tree(path)
        # An unobserved zero-height tip must still fail scaled validation.
        for metric in (pse, psc)
            @test_throws ArgumentError metric([1 2 0], tree)
            @test isfinite(only(getproperty(metric([1 2 0], tree; scale_vcv = false), nameof(metric))))
        end
    end
end

@testset "PSD shares sparse site trees and preserves the complete variance pool" begin
    tree = species_comb(20_000)
    n = 8
    matrix = sparse(repeat(1:n, inner = 3),
        [1 + mod(197i, tree.ntips) for i in 1:(3n)], repeat([1, 2, 3], n), n, tree.ntips)
    matrix[1, :] .= 0
    matrix[2, :] .= 0
    matrix[2, end] = 1
    matrix[3, findall(>(0), matrix[3, :])] .= 1
    for scaled in (false, true), compute_var in (false, true)
        options = (; scale_vcv = scaled, compute_var)
        combined = psd(matrix, tree; options...)
        variability = psv(matrix, tree; options...)
        richness = variability.richness
        expected = (variability.psv, variability.variance, psc(matrix, tree; scale_vcv = scaled).psc,
            variability.psv .* richness, variability.variance .* richness .^ 2,
            pse(matrix, tree; scale_vcv = scaled).pse, richness)
        for (field, value) in zip(propertynames(combined), expected)
            @test species_match(getproperty(combined, field), value)
        end
        @test all(getproperty(combined, a) !== getproperty(combined, b)
            for a in propertynames(combined) for b in propertynames(combined) if a != b)
        @test @allocated(psd(matrix, tree; options...)) < 12_000_000
        separate = @allocated begin
            psv(matrix, tree; options...)
            psr(matrix, tree; options...)
            pse(matrix, tree; scale_vcv = scaled)
            psc(matrix, tree; scale_vcv = scaled)
        end
        @test @allocated(psd(matrix, tree; options...)) < separate
    end
end

@testset "Species site-local paths and storage on a deep comb" begin
    tree = species_comb(20_000)
    sites = 64
    matrix = sparse(repeat(1:sites, inner = 3),
        [1 + mod(313i, tree.ntips) for i in 1:(3sites)], repeat([1, 2, 3], sites), sites, tree.ntips)
    saved = copy(matrix)
    for scaled in (false, true)
        rows, plan = species_prepared(matrix, tree, scaled)
        @test plan.tree.nnodes <= 2nnz(matrix) - 1
        @test plan.index !== nothing
        workspace = EvoKits._SpeciesSiteWorkspace(plan.tree.nnodes)
        for tips in plan.rows
            local_tree, _, local_tips = EvoKits._species_site!(workspace, plan, tips)
            @test local_tree.nnodes <= 2length(tips) - 1
            @test length(local_tips) == length(tips)
        end
        full_rows, full_plan = species_prepared(matrix, tree, scaled; strategy = :full)
        e = pse(matrix, tree; scale_vcv = scaled).pse
        c = psc(matrix, tree; scale_vcv = scaled).psc
        @test species_match(e, EvoKits._species_pse_values(full_rows, full_plan))
        @test species_match(c, EvoKits._species_psc_values(full_plan))
        @test @allocated(pse(matrix, tree; scale_vcv = scaled)) < 12_000_000
        @test @allocated(psc(matrix, tree; scale_vcv = scaled)) < 12_000_000
    end
    @test matrix == saved
    # The union covers the full pool, while each individual query is sparse.
    tree = species_comb(1024)
    matrix = sparse(repeat(1:512, inner = 2), 1:1024, repeat([1, 3], 512), 512, 1024)
    for scaled in (false, true)
        rows, plan = species_prepared(matrix, tree, scaled)
        @test plan.tree.nnodes == tree.nnodes
        @test plan.index !== nothing
        for (site, tips) in enumerate(plan.rows)
            # Every pair here consists of consecutive original tips. The LCA
            # and path lengths are known directly from the comb construction.
            a, b = 2site - 1, 2site
            raw = (b == tree.ntips ? 2 : b - a + 2) / 2
            expected_e = scaled ? (1 - (a - 1) / sqrt(a * (b == tree.ntips ? b - 1 : b))) * 0.75 : raw * 0.75
            expected_c = scaled ? expected_e / 0.75 : 2.0 - a
            local_tree, geometry, local_tips = EvoKits._species_site!(EvoKits._SpeciesSiteWorkspace(plan.tree.nnodes), plan, tips)
            @test local_tree.nnodes == 3
            @test EvoKits._pse_moments!(local_tree, local_tips, rows.weights[site], geometry,
                zeros(3), zeros(3)) ≈ expected_e
            @test EvoKits._psc_value!(EvoKits._ComdistntWorkspace(3), local_tree, local_tips, geometry) ≈ expected_c
        end
    end
end
