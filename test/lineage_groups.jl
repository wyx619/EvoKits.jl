using SparseArrays
using CSV

function lineage_test_tree(newick)
    return mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, newick)
        read_tree(path)
    end
end

# Independent root-depth arithmetic in BigFloat, followed by separate ancestor
# walks from each tip. No longest-descendant ages, preorder owners or event
# degree counting from the implementation are used.
function lineage_reference_ages(tree)
    incoming = Dict(Int(tree.child_of_edge[e]) => BigFloat(tree.edge_length[e]) for e in 1:tree.nedges)
    depths = zeros(BigFloat, tree.nnodes)
    for node in 1:tree.nnodes
        current = node
        while current != tree.root
            depths[node] += incoming[current]
            current = Int(tree.parent_of_node[current])
        end
    end
    return maximum(depths[tree.tip_ids]) .- depths
end

function lineage_reference_partition(tree, ages, cut, boundary)
    owners = Int[]
    for tip in tree.tip_ids
        owner = Int(tip)
        parent = Int(tree.parent_of_node[owner])
        while parent != 0 && (ages[parent] < cut || (boundary === :older && ages[parent] == cut))
            owner = parent
            parent = Int(tree.parent_of_node[parent])
        end
        push!(owners, owner)
    end
    groups = unique(owners)
    return [findfirst(==(owner), groups) for owner in owners], groups
end

function lineage_reference_target(tree, ages, k)
    if k == 1
        return lineage_reference_partition(tree, ages, ages[Int(tree.root)], :older)
    end
    # Enumerate candidate times, computing each partition from scratch. Unary
    # nodes are harmless candidates; no precomputed branching increments.
    events = sort!(unique(ages[.!tree.is_tip]); rev = true)
    for age in events
        membership, nodes = lineage_reference_partition(tree, ages, age, :younger)
        maximum(membership) >= k && return membership, nodes
    end
    error("no target partition")
end

function lineage_random_newick(rng, ntips)
    function emit(lo, hi, age, parent_age)
        if lo == hi
            return "t$lo:$(parent_age)"
        end
        mid = rand(rng, lo:(hi - 1))
        left_age = lo == mid ? 0 : rand(rng, 0:age)
        right_age = mid + 1 == hi ? 0 : rand(rng, 0:age)
        return "(" * emit(lo, mid, left_age, age) * "," *
            emit(mid + 1, hi, right_age, age) * "):$(parent_age - age)"
    end
    return emit(1, ntips, ntips * 4, ntips * 4) * ";"
end

@testset "Tree lineage groups: hand calculations, exact boundaries and ties" begin
    tree = lineage_test_tree("(((A:1,B:1):1,C:2):2,(D:2,E:2):2);")
    @test get_clades(tree; cut = 1).membership == [1, 1, 2, 3, 4]
    @test get_clades(tree; cut = 1, boundary = :younger).membership == collect(1:5)
    @test get_clades(tree; cut = 2).membership == [1, 1, 1, 2, 2]
    @test get_clades(tree; cut = 2, boundary = :younger).membership == [1, 1, 2, 3, 4]
    @test get_clades(tree; cut = 4).membership == ones(Int, 5)
    @test get_clades(tree; cut = 4, boundary = :younger).membership == [1, 1, 1, 2, 2]
    @test get_clades(tree; cut = 0).membership == collect(1:5)
    @test get_clades(tree; cut = prevfloat(2.0)).membership == [1, 1, 2, 3, 4]
    @test get_clades(tree; cut = nextfloat(2.0)).membership == [1, 1, 1, 2, 2]
    result = get_clades(tree; k = 3)
    @test result isa CladeGroupingResult
    @test result.requested_k == 3 && result.actual_k == 4
    @test result.cut == 2 && result.boundary === :younger && result.tree_height == 4
    @test result.membership == [1, 1, 2, 3, 4]
    @test result.group_labels == ["clade_1", "clade_2", "clade_3", "clade_4"]
    @test result.group_sizes == [2, 1, 1, 1]
    @test result.tip_labels == ["A", "B", "C", "D", "E"]
    for k in 1:5
        groups = get_clades(tree; k)
        @test get_clades(tree; cut = groups.cut, boundary = groups.boundary).membership == groups.membership
        @test groups.actual_k >= k && sum(groups.group_sizes) == 5
        @test groups.actual_k == length(groups.group_nodes) == length(groups.group_labels)
    end
    star = lineage_test_tree("(A:2,B:2,C:2,D:2);")
    @test get_clades(star; k = 2).actual_k == 4
    @test get_clades(star; k = 3).actual_k == 4
    @test get_clades(star; cut = 2).actual_k == 1
    zero = lineage_test_tree("((A:0,B:0):0,C:0);")
    @test get_clades(zero; cut = 0).actual_k == 1
    @test get_clades(zero; cut = 0, boundary = :younger).actual_k == 3
    @test get_clades(zero; k = 1).actual_k == 1
    @test get_clades(zero; k = 2).actual_k == 3
    terminal = lineage_test_tree("((A:0,B:0):2,C:2);")
    @test get_clades(terminal; cut = 0).membership == [1, 1, 2]
    @test get_clades(terminal; cut = 0, boundary = :younger).membership == [1, 2, 3]
    @test get_clades(terminal; k = 3).membership == [1, 2, 3]
    unary = lineage_test_tree("(((A:1,B:1):0):1,C:2);")
    @test get_clades(unary; k = 2).membership == [1, 1, 2]
    @test get_clades(unary; k = 3).membership == [1, 2, 3]
    singleton = lineage_test_tree("(A:0);")
    @test get_clades(singleton; cut = 0).membership == [1]
    @test get_clades(singleton; k = 1).group_sizes == [1]
    @test get_clades(singleton; cut = 0, boundary = :younger).actual_k == 1
    @test get_clades(tree; cut = 2).membership == get_clades(tree; cut = 2, atol = 1.0).membership
end

@testset "Lineage grouping and communities: saved local R references" begin
    fixture(name) = CSV.read(joinpath(@__DIR__, "fixtures", "lineages", name), DataFrame; delim = '\t')
    trees = fixture("trees.tsv")
    partitions = fixture("partitions_r.tsv")
    communities = fixture("communities_r.tsv")
    for entry in eachrow(trees)
        tree = lineage_test_tree(entry.newick)
        reference = partitions[partitions.case .== entry.case, :]
        for request in unique(reference.request)
            rows = reference[reference.request .== request, :]
            grouping = first(rows.mode) == "cut" ? get_clades(tree; cut = first(rows.value)) :
                get_clades(tree; k = Int(first(rows.value)))
            @test grouping.tip_labels == rows.species
            @test grouping.membership == rows.group
            @test grouping.actual_k == maximum(rows.group)
            saved = communities[(communities.case .== entry.case) .& (communities.request .== request), :]
            isempty(saved) && continue
            input = reshape(saved.input, 4, tree.ntips)
            for statistic in (:presence, :richness, :abundance), storage in (input, sparse(input))
                result = collapse_range(storage, grouping; statistic)
                values = [result.community.matrix[saved.site[i], saved.group[i]] for i in 1:nrow(saved)]
                @test values == saved[!, statistic]
                @test result.community.species_labels == grouping.group_labels
                @test result.grouping.membership == rows.group
            end
        end
    end
end

@testset "Tree lineage groups: independent random ancestor partitions" begin
    rng = MersenneTwister(20261106)
    for trial in 1:45
        n = rand(rng, 2:24)
        tree = lineage_test_tree(lineage_random_newick(rng, n))
        ages = lineage_reference_ages(tree)
        cuts = unique(vcat([0.0, Float64(maximum(ages)), rand(rng) * Float64(maximum(ages))],
            Float64.(ages[.!tree.is_tip])))
        for cut in cuts, boundary in (:older, :younger)
            actual = get_clades(tree; cut, boundary)
            expected, nodes = lineage_reference_partition(tree, ages, BigFloat(cut), boundary)
            @test actual.membership == expected
            @test actual.group_nodes == nodes
            @test actual.group_sizes == [count(==(i), expected) for i in 1:maximum(expected)]
            @test actual.actual_k == maximum(expected) && actual.requested_k === nothing
            for (group, node) in enumerate(nodes)
                @test Set(tree.tip_labels[expected .== group]) == Set(descendant_tips(tree, node; labels = true))
            end
        end
        for k in 1:n
            actual = get_clades(tree; k)
            expected, nodes = lineage_reference_target(tree, ages, k)
            @test actual.membership == expected && actual.group_nodes == nodes
            @test actual.actual_k >= k && actual.requested_k == k
            @test get_clades(tree; cut = actual.cut, boundary = actual.boundary).membership == expected
        end
    end
end

@testset "Tree lineage groups: validation, scaling and snapshot ownership" begin
    tree = lineage_test_tree("((A:1,B:1):2,C:3);")
    for cut in (-1, NaN, Inf, true, "one", big"1e400", big"1e-400", 3.1)
        @test_throws ArgumentError get_clades(tree; cut)
    end
    for k in (-1, 0, 4, 1.5, true, "two")
        @test_throws ArgumentError get_clades(tree; k)
    end
    @test_throws ArgumentError get_clades(tree)
    @test_throws ArgumentError get_clades(tree; cut = 1, k = 2)
    @test_throws ArgumentError get_clades(tree; cut = 1, boundary = :wrong)
    @test_throws ArgumentError get_clades(tree; k = 2, boundary = :older)
    for tolerance in (-1, Inf, NaN, true)
        @test_throws ArgumentError get_clades(tree; cut = 1, atol = tolerance)
        @test_throws ArgumentError get_clades(tree; cut = 1, rtol = tolerance)
    end
    @test_throws ArgumentError get_clades(lineage_test_tree("(A:1,B:2);" ); cut = 1)
    near = lineage_test_tree("(A:1,B:1.000000001);")
    @test get_clades(near; cut = 0).actual_k == 2
    @test_throws ArgumentError get_clades(near; cut = 0, rtol = 0)
    @test get_clades(near; cut = 0, rtol = 0, atol = 1e-8).actual_k == 2
    for bad_length in (-1.0, Inf, NaN)
        bad = deepcopy(tree)
        bad.edge_length[1] = bad_length
        @test_throws ArgumentError get_clades(bad; cut = 1)
    end
    overflow = lineage_test_tree("((A:1e308,B:1e308):1e308,C:1e308);")
    @test_throws ArgumentError get_clades(overflow; cut = 0)
    for factor in (1e-300, 1e300)
        scaled = deepcopy(tree)
        scaled.edge_length .*= factor
        # The cached root depths are deliberately not updated.
        @test get_clades(scaled; cut = factor).membership == get_clades(tree; cut = 1).membership
        @test get_clades(scaled; k = 2).membership == get_clades(tree; k = 2).membership
    end
    # Short young edges must survive a long root stem.
    stem = lineage_test_tree("((A:1,B:1):1e20);")
    @test get_clades(stem; cut = 1).actual_k == 1
    @test get_clades(stem; cut = prevfloat(1.0)).actual_k == 2
    result = get_clades(tree; cut = 1)
    saved = deepcopy(tree)
    result.tip_labels[1] = "changed"
    result.membership .= 99
    result.group_nodes .= 99
    result.group_sizes .= 99
    result.group_labels[1] = "changed"
    @test tree.tip_labels == saved.tip_labels && tree.edge_length == saved.edge_length
    @test get_clades(tree; cut = 1).membership == [1, 1, 2]
    labels = get_clades(tree; cut = 1).tip_labels
    tree.tip_labels[1] = "renamed"
    @test labels == ["A", "B", "C"]
    tree.tip_labels[1] = "B"
    @test_throws ArgumentError get_clades(tree; cut = 1)
end

function lineage_reference_community(matrix, species, grouping, statistic, normalize)
    result = zeros(BigFloat, size(matrix, 1), grouping.actual_k)
    for j in axes(matrix, 2), i in axes(matrix, 1)
        matrix[i, j] > 0 || continue
        tip = findfirst(==(species[j]), grouping.tip_labels)
        group = grouping.membership[tip]
        if statistic === :presence
            result[i, group] = 1
        else
            result[i, group] += statistic === :abundance ? BigFloat(matrix[i, j]) : 1
        end
    end
    if normalize
        for i in axes(result, 1)
            total = sum(result[i, :])
            total > 0 && (result[i, :] ./= total)
        end
    end
    return Float64.(result)
end

@testset "Lineage community collapse: mapping, statistics and beta integration" begin
    tree = lineage_test_tree("(((A:1,B:1):1,C:2):2,(D:2,E:2):2);")
    grouping = get_clades(tree; cut = 2)
    matrix = [2 4 0 0 0; 0 0 5 1 2; 0 0 0 0 0; 0 0 0 3 0]
    samples = ["ab", "cde", "empty", "d"]
    for storage in (matrix, sparse(matrix))
        for statistic in (:presence, :richness, :abundance), normalize in (false, true)
            data = CommunityData(storage; species_labels = tree.tip_labels, sample_labels = samples)
            result = collapse_range(data, tree; cut = 2, statistic, normalize)
            @test result isa LineageCommunityResult
            @test result.community.matrix isa SparseMatrixCSC
            expected = lineage_reference_community(matrix, tree.tip_labels, grouping, statistic, normalize)
            @test same_values(Matrix(result.community.matrix), expected)
            @test result.community.sample_labels == samples
            @test result.community.species_labels == grouping.group_labels
            @test result.grouping.membership == grouping.membership
            @test result.statistic === statistic && result.normalize === normalize
            reused = collapse_range(data, grouping; statistic, normalize)
            direct = collapse_range(storage, tree; cut = 2, statistic, normalize, sample_labels = samples)
            mapped = collapse_range(storage, grouping; statistic, normalize, sample_labels = samples)
            @test reused.community.matrix == direct.community.matrix == mapped.community.matrix == result.community.matrix
            @test beta_diss(result.community).sample_labels == samples
            @test isempty(nonzeros(result.community.matrix[3:3, :]))
        end
    end
    @test collapse_range(matrix, tree; cut = 2).community.matrix == [1 0; 1 1; 0 0; 0 1]
    @test collapse_range(matrix, tree; cut = 2, statistic = :richness).community.matrix == [2 0; 1 2; 0 0; 0 1]
    @test collapse_range(matrix, tree; cut = 2, statistic = :abundance).community.matrix == [6 0; 5 3; 0 0; 0 3]
    subset = CommunityData(sparse(matrix[:, [3, 1]]); species_labels = ["C", "A"], sample_labels = samples)
    result = collapse_range(subset, grouping)
    @test size(result.community) == (4, 2)
    @test nnz(result.community.matrix[:, 2:2]) == 0
    @test result.community.matrix[:, 1] == [1, 1, 0, 0]
    @test result.grouping.group_sizes == [3, 2] # Full tree, not observed subset.
    @test collapse_range(subset.matrix, tree; cut = 2, tip_labels = subset.species_labels).community.matrix == result.community.matrix
    reordered = CommunityData(sparse(matrix[:, [5, 2, 4, 1, 3]]); species_labels = tree.tip_labels[[5, 2, 4, 1, 3]])
    @test collapse_range(reordered, grouping).community.matrix == collapse_range(matrix, grouping).community.matrix
    nonempty = CommunityData(sparse(matrix[[1, 2, 4], :]); species_labels = tree.tip_labels)
    regions = phyloregion(beta_diss(collapse_range(nonempty, grouping).community); k = 2)
    @test length(regions.membership) == 3 && sum(regions.cluster_sizes) == 3
    @test collapse_range(zeros(Int, 0, 5), grouping).community.matrix == spzeros(Int, 0, 2)
    absent = CommunityData(spzeros(Int, 3, 0); species_labels = String[])
    @test size(collapse_range(absent, grouping).community) == (3, 2)
    bool_result = collapse_range(matrix .> 0, grouping; statistic = :abundance)
    @test eltype(bool_result.community.matrix) === Int
end

@testset "Lineage collapse: independent randomized sparse aggregation" begin
    rng = MersenneTwister(261106)
    for trial in 1:35
        p, n = rand(rng, 2:20), rand(rng, 2:16)
        tree = lineage_test_tree(lineage_random_newick(rng, p))
        grouping = get_clades(tree; k = rand(rng, 1:p))
        matrix = rand(rng, 0:5, n, p)
        matrix[rand(rng, n, p) .> 0.2] .= 0
        columns = randperm(rng, p)[1:rand(rng, 1:p)]
        species = tree.tip_labels[columns]
        input = matrix[:, columns]
        saved = copy(input)
        for statistic in (:presence, :richness, :abundance), normalize in (false, true)
            data = CommunityData(sparse(input); species_labels = species)
            expected = lineage_reference_community(input, species, grouping, statistic, normalize)
            actual = collapse_range(data, grouping; statistic, normalize)
            dense = collapse_range(input, grouping; statistic, normalize, tip_labels = species)
            @test same_values(Matrix(actual.community.matrix), expected)
            @test actual.community.matrix == dense.community.matrix
            @test actual.community.species_labels == grouping.group_labels
            @test actual.grouping.group_nodes == grouping.group_nodes
            @test nnz(actual.community.matrix) <= count(>(0), input)
            if normalize
                @test all(value -> value == 0 || isapprox(value, 1), vec(sum(actual.community.matrix; dims = 2)))
            elseif statistic === :abundance
                @test sum(actual.community.matrix) == sum(input)
            elseif statistic === :richness
                @test sum(actual.community.matrix) == count(>(0), input)
            end
        end
        @test input == saved
    end
end

@testset "Lineage collapse: invalid inputs, overflow and result isolation" begin
    tree = lineage_test_tree("((A:1,B:1):1,C:2);")
    grouping = get_clades(tree; cut = 1)
    data = CommunityData(sparse([1 2 0; 0 0 0]); species_labels = tree.tip_labels)
    @test_throws ArgumentError collapse_range(data, grouping; statistic = :wrong)
    @test_throws ArgumentError collapse_range(zeros(2, 2), grouping)
    @test_throws ArgumentError collapse_range(zeros(2, 3), grouping; tip_labels = ["A", "B", "X"])
    @test_throws ArgumentError collapse_range(CommunityData(zeros(2, 1); species_labels = ["X"]), grouping)
    @test_throws ArgumentError collapse_range([-1 0 0], grouping)
    @test_throws ArgumentError collapse_range([NaN 0 0], grouping)
    @test_throws ArgumentError collapse_range([typemax(Int) 1 0], grouping; statistic = :abundance)
    @test_throws ArgumentError collapse_range([1e308 1e308 0.0], grouping; statistic = :abundance)
    @test eltype(collapse_range(Int32[2 4 0], grouping; statistic = :abundance).community.matrix) === Int32
    for mutation in (:membership, :size, :labels)
        damaged = deepcopy(grouping)
        mutation === :membership && (damaged.membership[1] = 99)
        mutation === :size && (damaged.group_sizes[1] = 99)
        mutation === :labels && (damaged.group_labels[1] = "")
        @test_throws ArgumentError collapse_range(data, damaged)
    end
    result = collapse_range(data, grouping)
    result.grouping.tip_labels[1] = "changed"
    result.grouping.membership .= 99
    result.grouping.group_labels[1] = "changed"
    result.grouping.group_nodes .= 99
    result.grouping.group_sizes .= 99
    result.community.sample_labels[1] = "changed"
    result.community.species_labels[1] = "changed"
    nonzeros(result.community.matrix) .= 99
    @test grouping.tip_labels == tree.tip_labels && grouping.membership == [1, 1, 2]
    @test grouping.group_labels == ["clade_1", "clade_2"] && grouping.group_sizes == [2, 1]
    @test data.sample_labels == ["sample_1", "sample_2"] && data.matrix == [1 2 0; 0 0 0]
    @test collapse_range(data, grouping).community.matrix == [1 0; 0 0]
end

function lineage_balanced_newick(p)
    function emit(lo, hi, age, parent_age)
        lo == hi && return "t$lo:$parent_age"
        mid = (lo + hi) ÷ 2
        return "(" * emit(lo, mid, age - 1, age) * "," *
            emit(mid + 1, hi, age - 1, age) * "):$(parent_age - age)"
    end
    height = ceil(Int, log2(p))
    return emit(1, p, height, height) * ";"
end

@testset "Lineage grouping and aggregation avoid dense or descendant-set storage" begin
    n, p = 2000, 8192
    tree = lineage_test_tree(lineage_balanced_newick(p))
    grouping = get_clades(tree; k = 50)
    @test grouping.actual_k == 64
    get_clades(tree; cut = grouping.cut)
    @test (@allocated get_clades(tree; cut = grouping.cut)) < 5_000_000
    @test Base.summarysize(grouping) < 1_000_000
    rng = MersenneTwister(604)
    matrix = sparse(rand(rng, 1:n, 6000), rand(rng, 1:p, 6000), ones(Int, 6000), n, p)
    data = CommunityData(matrix; species_labels = tree.tip_labels)
    collapse_range(data, grouping; statistic = :abundance)
    @test (@allocated collapse_range(data, grouping; statistic = :abundance)) < 25_000_000
    result = collapse_range(data, grouping; statistic = :abundance)
    @test result.community.matrix isa SparseMatrixCSC && nnz(result.community.matrix) <= nnz(matrix)
    @test sum(result.community.matrix) == sum(matrix)
    @test Base.summarysize(result) < 2_000_000
end
