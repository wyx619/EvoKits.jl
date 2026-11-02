using CSV
using LinearAlgebra

# Independent global-minimum UPGMA. Each candidate is recomputed from original
# site pairs in BigFloat; it does not use packed distances or neighbor chains.
function region_reference(distances)
    n = size(distances, 1)
    clusters = [[i] for i in 1:n]
    nodes = -collect(1:n)
    merge, heights = zeros(Int, n - 1, 2), zeros(n - 1)
    for step in 1:(n - 1)
        best, left, right = BigFloat(Inf), 0, 0
        for j in 2:length(clusters), i in 1:(j - 1)
            value = sum(BigFloat(distances[a, b]) for a in clusters[i], b in clusters[j]) /
                (length(clusters[i]) * length(clusters[j]))
            if value < best
                best, left, right = value, i, j
            end
        end
        merge[step, :] = [nodes[left], nodes[right]]
        heights[step] = Float64(best)
        append!(clusters[left], clusters[right])
        sort!(clusters[left])
        nodes[left] = step
        deleteat!(clusters, right)
        deleteat!(nodes, right)
    end
    return merge, heights
end

function region_reference_cut(merge, n, k)
    descendants = Dict(-i => [i] for i in 1:n)
    clusters = Dict(-i => [i] for i in 1:n)
    for step in 1:(n - k)
        a, b = merge[step, 1], merge[step, 2]
        descendants[step] = vcat(descendants[a], descendants[b])
        delete!(clusters, a)
        delete!(clusters, b)
        clusters[step] = descendants[step]
    end
    membership = zeros(Int, n)
    for (region, members) in enumerate(sort!(collect(values(clusters)); by = minimum))
        membership[members] .= region
    end
    return membership
end

function region_reference_summary(distances, membership)
    k = maximum(membership)
    region_distances, ed = zeros(k, k), fill(NaN, k)
    for a in 1:k, b in (a + 1):k
        first, second = findall(==(a), membership), findall(==(b), membership)
        value = sum(BigFloat(distances[i, j]) for i in first, j in second) / (length(first) * length(second))
        region_distances[a, b] = region_distances[b, a] = Float64(value)
    end
    if k > 1
        for a in 1:k
            ed[a] = Float64(sum(BigFloat(region_distances[a, b]) for b in 1:k if a != b) / (k - 1))
        end
    end
    return region_distances, ed
end

# Verify an entire hierarchy is a valid global-minimum UPGMA hierarchy even
# when exact ties admit several different, equally valid trees.
function check_region_hierarchy(distances, result)
    n = size(distances, 1)
    active = Dict(-i => [i] for i in 1:n)
    @test size(result.merge) == (n - 1, 2)
    @test length(result.heights) == n - 1 && issorted(result.heights)
    for step in 1:(n - 1)
        a, b = result.merge[step, 1], result.merge[step, 2]
        @test a != b && a < step && b < step && haskey(active, a) && haskey(active, b)
        best = minimum(sum(BigFloat(distances[i, j]) for i in left, j in right) /
            (length(left) * length(right)) for (x, left) in active for (y, right) in active if x < y)
        selected = sum(BigFloat(distances[i, j]) for i in active[a], j in active[b]) /
            (length(active[a]) * length(active[b]))
        @test result.heights[step] ≈ Float64(best) rtol = 1e-12 atol = 1e-14
        @test result.heights[step] ≈ Float64(selected) rtol = 1e-12 atol = 1e-14
        active[step] = vcat(active[a], active[b])
        delete!(active, a)
        delete!(active, b)
    end
    @test length(active) == 1 && sort!(only(values(active))) == collect(1:n)
end

@testset "Regionalization hand calculations, boundaries and snapshots" begin
    distances = [0.0 1 8 10; 1 0 6 12; 8 6 0 2; 10 12 2 0]
    labels = ["AB_1", "AB_2", "CD_1", "CD_2"]
    saved = copy(distances)
    result = phyloregion(distances; k = 2, sample_labels = labels)
    @test result isa PhyloregionResult
    @test result.membership == [1, 1, 2, 2]
    @test result.cluster_sizes == [2, 2]
    @test result.region_distances == [0 9; 9 0]
    @test result.distinctiveness == [9, 9]
    @test result.merge == [-1 -2; -3 -4; 1 2]
    @test result.heights == [1, 2, 9]
    @test result.k == 2 && result.method === :average && result.component === :distance
    @test result.sample_labels == labels
    result3 = phyloregion(distances; k = 3)
    @test result3.membership == [1, 1, 2, 3]
    @test result3.region_distances == [0 7 11; 7 0 2; 11 2 0]
    # Equal weight per other region, rather than per site in those regions.
    @test result3.distinctiveness == [9, 4.5, 6.5]
    @test phyloregion(distances; k = 1).membership == ones(Int, 4)
    @test isnan(only(phyloregion(distances; k = 1).distinctiveness))
    @test phyloregion(distances; k = 4).membership == collect(1:4)
    @test phyloregion(distances; k = 4).region_distances == distances
    singleton = phyloregion(zeros(1, 1); k = 1, sample_labels = [123])
    @test singleton.membership == singleton.cluster_sizes == [1]
    @test singleton.sample_labels == ["123"]
    @test size(singleton.merge) == (0, 2) && isempty(singleton.heights)
    @test singleton.region_distances == zeros(1, 1) && isnan(only(singleton.distinctiveness))
    @test phyloregion(sparse(distances); k = 2).membership == result.membership
    @test phyloregion(view(distances, :, :); k = 2).merge == result.merge
    # The input remains independent of all result arrays, including labels.
    result.sample_labels[1] = "changed"
    result.membership[1] = 99
    result.cluster_sizes[1] = 99
    result.region_distances[1, 2] = 99
    result.distinctiveness[1] = 99
    result.heights[1] = 99
    result.merge[1, 1] = 99
    @test distances == saved && labels == ["AB_1", "AB_2", "CD_1", "CD_2"]
    @test phyloregion(distances; k = 2).heights == [1, 2, 9]
    check_region_hierarchy(distances, phyloregion(distances; k = 2))
end

@testset "Regionalization defined ties and independent random hierarchies" begin
    equal = ones(5, 5) - Matrix{Float64}(I, 5, 5)
    first = phyloregion(equal; k = 2)
    @test first.merge == [-1 -2; 1 -3; 2 -4; 3 -5]
    @test first.membership == [1, 1, 1, 1, 2]
    @test first.heights == ones(4)
    @test phyloregion(zeros(5, 5); k = 3).membership == [1, 1, 1, 2, 3]
    rng = MersenneTwister(20261103)
    for trial in 1:40
        n = rand(rng, 2:16)
        # Half the cases have exact ties. Do not impose the slow reference's
        # different global tie rule on nearest-neighbor-chain clustering.
        distances = rand(rng, n, n)
        trial <= 20 && (distances .= rand(rng, 0:5, n, n))
        distances = Matrix(Symmetric(distances))
        distances[diagind(distances)] .= 0
        result = phyloregion(distances; k = rand(rng, 1:n))
        check_region_hierarchy(distances, result)
        @test result.merge == phyloregion(distances; k = 1).merge
        for k in unique([1, n, result.k])
            actual = phyloregion(distances; k)
            expected = region_reference_cut(actual.merge, n, k)
            regions, ed = region_reference_summary(distances, expected)
            @test actual.membership == expected
            @test actual.cluster_sizes == [count(==(i), expected) for i in 1:k]
            @test same_values(actual.region_distances, regions)
            @test same_values(actual.distinctiveness, ed)
        end
        if trial > 20
            reference, heights = region_reference(distances)
            @test same_values(result.heights, heights)
            for k in 1:n
                @test phyloregion(distances; k).membership == region_reference_cut(reference, n, k)
            end
            permutation = randperm(rng, n)
            reordered = phyloregion(distances[permutation, permutation]; k = result.k)
            @test (reordered.membership .== reordered.membership') ==
                (result.membership[permutation] .== result.membership[permutation]')
        end
    end
end

@testset "Regionalization validation and extreme finite distances" begin
    distances = [0.0 1 3; 1 0 2; 3 2 0]
    for k in (0, -1, 4, 1.5, true)
        @test_throws ArgumentError phyloregion(distances; k)
    end
    @test_throws ArgumentError phyloregion(distances; k = 2, method = :ward)
    for bad in (zeros(0, 0), zeros(2, 3), [0 1; 2 0], [1 2; 2 0],
        [0 -1; -1 0], [0 NaN; NaN 0], [0 Inf; Inf 0], [NaN 1; 1 0])
        @test_throws ArgumentError phyloregion(bad; k = 1)
    end
    for labels in (["A", "A", "B"], ["", "B", "C"], ["A", "B"], [1, missing, 3])
        @test_throws ArgumentError phyloregion(distances; k = 2, sample_labels = labels)
    end
    for value in (big"1e400", big"1e-400")
        bad = [BigFloat(0) value; value BigFloat(0)]
        @test_throws ArgumentError phyloregion(bad; k = 1)
    end
    for factor in (1e-300, 1e300)
        result = phyloregion(distances .* factor; k = 2)
        @test result.membership == [1, 1, 2]
        @test result.heights ./ factor ≈ [1, 2.5]
        @test result.distinctiveness ./ factor ≈ [2.5, 2.5]
    end
    huge = fill(1e308, 4, 4)
    huge[diagind(huge)] .= 0
    result = phyloregion(huge; k = 3)
    @test all(==(1e308), result.heights)
    @test all(==(1e308), result.distinctiveness)
    @test result.region_distances[1, 2] == 1e308
    # Sorting discoveries must put independent low merges before earlier-found
    # high ones, and translate parent references to the sorted hierarchy.
    out_of_order = [0 9 30 30; 9 0 30 30; 30 30 0 1; 30 30 1 0]
    result = phyloregion(out_of_order; k = 3)
    @test result.heights == [1, 9, 30]
    @test result.membership == [1, 2, 3, 3]
    check_region_hierarchy(out_of_order, result)
end

@testset "Regionalization integrates full taxonomic and phylogenetic beta" begin
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "((A:1,B:2):3,(C:4,D:2):1);")
        tree = read_tree(path)
        data = CommunityData(sparse([1 1 0 0; 1 0 1 0; 0 1 1 1; 0 0 1 1]);
            species_labels = tree.tip_labels, sample_labels = ["AB", "AC", "BCD", "CD"])
        for beta in (beta_diss(data), phylobeta(data, tree)), component in (:total, :turnover, :nestedness)
            result = phyloregion(beta; k = 2, component)
            matrix = phyloregion(getproperty(beta, component); k = 2, sample_labels = beta.sample_labels)
            @test result.component === component
            @test result.sample_labels == data.sample_labels
            @test result.membership == matrix.membership
            @test result.merge == matrix.merge && result.heights == matrix.heights
            @test result.region_distances == matrix.region_distances
            @test result.distinctiveness == matrix.distinctiveness
        end
        @test_throws ArgumentError phyloregion(beta_diss(data); k = 2, component = :wrong)
        @test_throws ArgumentError phyloregion(beta_diss(data; pairs = [(1, 2)]); k = 1)
        @test_throws ArgumentError phyloregion(phylobeta(data, tree; pairs = [(1, 2)]); k = 1)
        with_empty = vcat(data.matrix, spzeros(Int, 1, 4))
        @test_throws ArgumentError phyloregion(beta_diss(with_empty); k = 2)
        @test_throws ArgumentError phyloregion(phylobeta(with_empty, tree); k = 2)
    end
end

@testset "Regionalization saved local R hierarchy and summaries for every k" begin
    inputs = CSV.read(joinpath(@__DIR__, "fixtures", "regionalization_input.tsv"), DataFrame; delim = '\t')
    r_merges = CSV.read(joinpath(@__DIR__, "fixtures", "regionalization_merges_r.tsv"), DataFrame; delim = '\t')
    r_members = CSV.read(joinpath(@__DIR__, "fixtures", "regionalization_membership_r.tsv"), DataFrame; delim = '\t')
    r_regions = CSV.read(joinpath(@__DIR__, "fixtures", "regionalization_regions_r.tsv"), DataFrame;
        delim = '\t', missingstring = "NaN")
    for case in unique(inputs.case)
        input = inputs[inputs.case .== case, :]
        merges = r_merges[r_merges.case .== case, :]
        n = maximum(input.i)
        distances = reshape(input.distance, n, n)
        labels = unique(input.label_i)
        reference_merge = hcat(merges.left, merges.right)
        for k in 1:n
            members = r_members[(r_members.case .== case) .& (r_members.k .== k), :]
            regions = r_regions[(r_regions.case .== case) .& (r_regions.k .== k), :]
            result = phyloregion(distances; k, sample_labels = labels)
            @test result.sample_labels == members.label
            @test result.membership == members.membership
            @test result.membership == region_reference_cut(reference_merge, n, k)
            @test result.cluster_sizes == [count(==(i), members.membership) for i in 1:k]
            @test same_values(result.heights, merges.height)
            @test same_values(result.region_distances, reshape(regions.distance, k, k))
            @test same_values(result.distinctiveness, Float64.(coalesce.(regions.distinctiveness[1:k], NaN)))
        end
    end
end

@testset "Regionalization packed workspace avoids a full distance copy" begin
    rng = MersenneTwister(777)
    n = 1000
    distances = Matrix(Symmetric(rand(rng, n, n)))
    distances[diagind(distances)] .= 0
    phyloregion(distances; k = 8)
    bytes = @allocated phyloregion(distances; k = 8)
    # The caller already holds 8 MB; one packed work triangle is about 4 MB.
    @test bytes < 7_000_000
    result = phyloregion(distances; k = 8)
    @test sum(result.cluster_sizes) == n && length(result.cluster_sizes) == 8
    @test Base.summarysize(result) < 200_000
end
