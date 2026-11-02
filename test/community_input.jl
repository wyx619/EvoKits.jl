@testset "CommunityData matrix snapshots, labels and storage" begin
    matrix = [0.0 1e-20 3; 2 0 0]
    samples, species = [11, 22], ["C", "A", "B"]
    for input in (matrix, sparse(matrix), matrix .> 0, view(matrix, :, :))
        saved = copy(input)
        data = CommunityData(input; sample_labels = samples, species_labels = species)
        @test data isa CommunityData
        @test size(data) == size(data.matrix) == (2, 3)
        @test size(data, 1) == 2 && size(data, 3) == 1
        @test data.sample_labels == ["11", "22"]
        @test data.species_labels == species
        @test data.species_labels !== species
        @test eltype(data.matrix) === eltype(input)
        @test issparse(data.matrix) == issparse(input)
        @test data.matrix == saved
        data.matrix[1, 1] = 1
        data.species_labels[1] = "changed"
        @test input == saved && species == ["C", "A", "B"]
        forced = CommunityData(input; species_labels = species, storage = :sparse)
        @test forced.matrix isa SparseMatrixCSC
        @test forced.matrix == input
    end
    explicit = sparse([1, 2, 2], [1, 1, 2], [0.0, nextfloat(0.0), 2.0], 2, 3)
    saved = copy(explicit)
    data = CommunityData(explicit; storage = :sparse)
    @test nnz(data.matrix) == 2 && nnz(explicit) == 3
    @test nonzeros(explicit) == nonzeros(saved)
    nonzeros(explicit)[2] = 10
    @test data.matrix[2, 1] == nextfloat(0.0)
    @test CommunityData(zeros(0, 0)).sample_labels == String[]
    @test CommunityData(zeros(2, 0)).sample_labels == ["sample_1", "sample_2"]
    @test CommunityData(zeros(0, 2)).species_labels == ["species_1", "species_2"]
    @test CommunityData(ones(1, 1); species_labels = [SubString("taxon", 1, 3)]).species_labels == ["tax"]
    @test_throws ArgumentError CommunityData(matrix; storage = :wrong)
    for labels in (["A"], ["A", "A", "B"], ["A", " ", "B"], ["A", missing, "B"], [1, "1", 2])
        @test_throws ArgumentError CommunityData(matrix; species_labels = labels)
    end
    for labels in (["x"], ["x", "x"], ["", "y"], [true, false], [1.0, 2.0])
        @test_throws ArgumentError CommunityData(matrix; sample_labels = labels)
    end
    for bad in (-1.0, Inf, -Inf, NaN), input in ([bad 0], sparse([bad 0]))
        @test_throws ArgumentError CommunityData(input)
    end
    @test_throws ArgumentError CommunityData([1 missing])
    @test_throws ArgumentError CommunityData([1 "bad"])
end

@testset "CommunityData long tables and indexed triplets" begin
    table = DataFrame(site = [20, 10, 20, 30, 10], taxon = ["B", "A", "B", "C", "B"], abundance = [2, 3, 4, 0, 1])
    saved = deepcopy(table)
    data = CommunityData(table; sample = :site, species = :taxon, value = :abundance)
    @test data.sample_labels == ["20", "10", "30"]
    @test data.species_labels == ["B", "A", "C"]
    @test data.matrix isa SparseMatrixCSC{Int}
    @test data.matrix == [6 0 0; 1 3 0; 0 0 0]
    @test nnz(data.matrix) == 3
    @test isequal(table, saved)
    reordered = CommunityData(table; sample = "site", species = "taxon", value = "abundance",
        sample_labels = [10, 30, 20, 40], species_labels = ["A", "B", "C", "D"])
    expected = [3 1 0 0; 0 0 0 0; 0 6 0 0; 0 0 0 0]
    @test reordered.matrix == expected
    @test reordered.sample_labels == ["10", "30", "20", "40"]
    triplets = CommunityData([3, 1, 3, 2, 1], [2, 1, 2, 3, 2], table.abundance;
        sample_labels = [10, 30, 20, 40], species_labels = ["A", "B", "C", "D"])
    @test triplets.matrix == expected
    @test triplets.species_labels == reordered.species_labels
    occurrences = CommunityData(select(table, :site, :taxon); sample = :site, species = :taxon, value = nothing)
    @test occurrences.matrix == [2 0 0; 1 1 0; 0 0 1]
    bool_counts = CommunityData([1, 1], [1, 1], [true, true]; sample_labels = [1], species_labels = ["A"])
    @test bool_counts.matrix[1, 1] == 2
    @test eltype(bool_counts.matrix) === Int
    # Mixed column types are promoted after validating each record.
    mixed = CommunityData(DataFrame(sample = ["a", "b"], species = ["A", "B"], value = Any[2, 0.5]))
    @test mixed.matrix == [2 0; 0 0.5]
    empty = DataFrame(sample = String[], species = String[], value = Float64[])
    @test size(CommunityData(empty)) == (0, 0)
    @test size(CommunityData(empty; sample_labels = ["empty"], species_labels = ["A", "B"])) == (1, 2)
    @test size(CommunityData(Int[], Int[], Float64[]; sample_labels = String[], species_labels = ["A"])) == (0, 1)
    for policy in (:error, :wrong)
        @test_throws ArgumentError CommunityData(table; sample = :site, species = :taxon, value = :abundance, duplicates = policy)
    end
    @test_throws ArgumentError CommunityData(DataFrame(sample = ["a", "a"], species = ["A", "A"], value = [0, 0]); duplicates = :error)
    for kw in ((sample = :absent,), (species = :absent,), (value = :absent,),
        (sample_labels = [10, 20],), (species_labels = ["A", "B"],))
        @test_throws ArgumentError CommunityData(table; sample = :site, species = :taxon, value = :abundance, kw...)
    end
    for bad in (missing, -1, NaN, Inf, "2"), values in ([bad, 1], Any[bad, 1])
        @test_throws ArgumentError CommunityData([1, 1], [1, 1], values; sample_labels = ["a"], species_labels = ["A"])
    end
    for bad in (missing, "", " ", true, 1.5)
        @test_throws ArgumentError CommunityData(DataFrame(sample = [bad], species = ["A"], value = [1]))
        @test_throws ArgumentError CommunityData(DataFrame(sample = ["a"], species = [bad], value = [1]))
    end
    for (rows, cols, values) in (([0], [1], [1]), ([2], [1], [1]), ([1.0], [1], [1]),
        ([true], [1], [1]), ([1], [0], [1]), ([1], [2], [1]), ([1], [1], [1, 2]))
        @test_throws ArgumentError CommunityData(rows, cols, values; sample_labels = ["a"], species_labels = ["A"])
    end
    for values in ([typemax(Int), 1], [typemax(UInt), UInt(1)], [1e308, 1e308])
        @test_throws ArgumentError CommunityData([1, 1], [1, 1], values; sample_labels = ["a"], species_labels = ["A"])
    end
    large_count = CommunityData([1, 1], [1, 1], BigInt[typemax(Int), 1]; sample_labels = ["a"], species_labels = ["A"])
    @test large_count.matrix[1, 1] == BigInt(typemax(Int)) + 1
end

@testset "Tree alignment preserves the full tree and sparse abundance" begin
    mktempdir() do dir
        path = joinpath(dir, "input.nwk")
        write_tree(path, "((A:1,B:2):3,(C:4,D:5):6);")
        tree = read_tree(path)
        source = [7 2; 0 3; 0 0]
        padded = [2 0 7 0; 3 0 0 0; 0 0 0 0]
        samples = ["first", "second", "empty"]
        tree_before = deepcopy(tree)
        for input in (source, sparse(source))
            data = CommunityData(input; sample_labels = samples, species_labels = ["C", "A"])
            aligned = align_community(data, tree)
            @test aligned.matrix == padded
            @test aligned.species_labels == tree.tip_labels
            @test aligned.sample_labels == samples
            @test issparse(aligned.matrix) == issparse(input)
            @test aligned.sample_labels !== data.sample_labels
            @test data.matrix == source && data.species_labels == ["C", "A"]
            @test_throws ArgumentError align_community(data, tree; missing_tips = :error)
            @test_throws ArgumentError align_community(data, tree; missing_tips = :wrong)
            strict = align_community(aligned, tree; missing_tips = :error)
            @test strict.matrix == padded && strict.matrix !== aligned.matrix
            aligned.matrix[1, 1] = 999
            @test strict.matrix == padded && data.matrix == source
            @test same_values(relative_pd(data, tree).rpd, relative_pd(padded, tree).rpd)
            @test same_values(phylo_endemism(data, tree).pe, phylo_endemism(padded, tree).pe)
            unknown = CommunityData(zeros(1, 1); species_labels = ["X"])
            @test_throws ArgumentError align_community(unknown, tree)
            @test_throws ArgumentError pd(unknown, tree)
        end
        @test tree.tip_labels == tree_before.tip_labels
        @test tree.edge_length == tree_before.edge_length
        @test tree.parent_of_node == tree_before.parent_of_node
        empty = CommunityData(zeros(2, 0); sample_labels = ["x", "y"])
        @test align_community(empty, tree).matrix == zeros(2, 4)
        @test pd(empty, tree).pd == [0, 0]
        @test all(isnan, relative_pd(empty, tree).rpd)
        zero_rows = CommunityData(spzeros(0, 1); species_labels = ["A"])
        @test size(align_community(zero_rows, tree)) == (0, 4)
        @test isempty(phylo_endemism(zero_rows, tree).pe)
    end
end

@testset "CommunityData metrics and prepared context match explicit matrix input" begin
    rng = MersenneTwister(20261105)
    mktempdir() do dir
        path = joinpath(dir, "input.nwk")
        write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
        tree = read_tree(path)
        padded = [2 0 5 0; 0 0 3 0; 1 0 0 0; 0 0 0 0]
        samples = ["ac", "c", "a", "empty"]
        for input in (padded[:, [3, 1]], sparse(padded[:, [3, 1]]))
            data = CommunityData(input; sample_labels = samples, species_labels = ["C", "A"])
            # Tree-free diagnostics retain the original species-column order.
            check_context_result(weighted_endemism(data), weighted_endemism(input; sample_labels = samples))
            for family in (:sorensen, :jaccard), pairs in (nothing, [(2, 1), (4, 4), (1, 1)])
                check_context_result(beta_diss(data; index_family = family, pairs = pairs),
                    beta_diss(input; index_family = family, pairs = pairs, sample_labels = samples))
                check_context_result(phylobeta(data, tree; index_family = family, pairs = pairs),
                    phylobeta(padded, tree; index_family = family, pairs = pairs, sample_labels = samples))
            end
            for f in (pd, phylo_endemism, relative_pd)
                check_context_result(f(data, tree), f(padded, tree; sample_labels = samples))
                @test_throws ArgumentError f(data, tree; tip_labels = ["wrong"])
                @test_throws ArgumentError f(data, tree; sample_labels = samples)
            end
            @test isequal(pd_dataframe(data, tree), pd_dataframe(padded, tree; sample_labels = samples))
            @test same_values(pd(data, tree; include_root = false).pd, pd(padded, tree; include_root = false).pd)
            @test same_values(phylo_endemism(data, tree; weighted = false).pe, phylo_endemism(padded, tree; weighted = false).pe)
            context = prepare_community(data, tree)
            explicit = prepare_community(padded, tree; sample_labels = samples)
            check_context_result(context, explicit)
            for f in (pd, weighted_endemism, phylo_endemism, relative_pd, phylobeta)
                check_context_result(f(context), f(explicit))
            end
            for f in (unifrac, phylosor, comdist, comdistnt)
                @test same_values(f(data, tree), f(padded, tree))
            end
            @test same_values(comdist(data, tree; abundance_weighted = true), comdist(padded, tree; abundance_weighted = true))
            @test same_values(comdistnt(data, tree; abundance_weighted = true), comdistnt(padded, tree; abundance_weighted = true))
            for f in (raoD, psv, psr, pse, psc, psd, pd_calc)
                check_context_result(f(data, tree), f(padded, tree))
            end
            check_context_result(pd_bootstrap(data, tree; reps = 5, rng = MersenneTwister(33)),
                pd_bootstrap(padded, tree; reps = 5, rng = MersenneTwister(33)))
            for f in (ses_pd, ses_mpd, ses_mntd)
                check_context_result(f(data, tree; runs = 5, rng = MersenneTwister(12)),
                    f(padded, tree; runs = 5, rng = MersenneTwister(12), sample_labels = samples))
            end
            for f in (ses_pd_dataframe, ses_mpd_dataframe, ses_mntd_dataframe)
                @test isequal(f(data, tree; runs = 5, rng = MersenneTwister(12)),
                    f(padded, tree; runs = 5, rng = MersenneTwister(12), sample_labels = samples))
            end
            data.matrix[1, 1] = 100
            data.sample_labels[1] = "changed"
            check_context_result(context, explicit)
        end
        # Repeated calls on a deliberately aligned snapshot use its matrix
        # directly, while metric results still own their label/result arrays.
        aligned = CommunityData(sparse(padded); sample_labels = samples, species_labels = tree.tip_labels)
        @test EvoKits._community_tree_data(aligned, tree) === aligned
        result = pd(aligned, tree)
        result.sample_labels[1] = "changed"
        @test aligned.sample_labels == samples
        @test_throws ArgumentError beta_diss(aligned; sample_labels = samples)
        @test_throws ArgumentError weighted_endemism(aligned; tip_labels = tree.tip_labels)
        # Mutating internal arrays violates the contract, but is checked at the
        # next public entry instead of propagating invalid labels or values.
        aligned.species_labels[1] = "X"
        @test_throws ArgumentError pd(aligned, tree)
        @test_throws ArgumentError beta_diss(CommunityData(zeros(1, 1)); index_family = :wrong)
    end
end

@testset "Sparse input helpers and large CommunityData stay sparse" begin
    rng = MersenneTwister(20261106)
    matrix = rand(rng, 0:4, 12, 9)
    matrix[matrix .< 3] .= 0
    perm = randperm(rng, size(matrix, 2))
    for f in (EvoKits._tools_abundance_only, EvoKits._tools_presence_only, EvoKits._tools_presence_matrix)
        @test f(sparse(matrix), perm) == f(matrix, perm)
    end
    mktempdir() do dir
        path = joinpath(dir, "sparse.nwk")
        write_tree(path, "(A:1,B:1,C:1,D:1,E:1,F:1,G:1,H:1,I:1);")
        tree = read_tree(path)
        @test EvoKits._community_presence(sparse(matrix), tree, tree.tip_labels[invperm(perm)]) ==
            EvoKits._community_presence(matrix, tree, tree.tip_labels[invperm(perm)])
    end
    n, p = 1500, 12000
    rows, cols = rand(rng, 1:n, 4000), rand(rng, 1:p, 4000)
    samples, species = ["site_$i" for i in 1:n], ["s$i" for i in 1:p]
    build() = CommunityData(rows, cols, ones(Int, length(rows)); sample_labels = samples, species_labels = species)
    data = build()
    @test @allocated(build()) < 8_000_000
    @test Base.summarysize(data) < 4_000_000
    pairs = [(1, 1500), (3, 4)]
    beta_diss(data; pairs = pairs)
    @test @allocated(beta_diss(data; pairs = pairs)) < 8_000_000
    @test nnz(data.matrix) <= 4000
end

@testset "CommunityData randomized long records, permutations and zero padding" begin
    rng = MersenneTwister(20261107)
    mktempdir() do dir
        path = joinpath(dir, "random.nwk")
        for trial in 1:30
            n, p = rand(rng, 2:12), rand(rng, 3:25)
            write_tree(path, random_community_tree(rng, p))
            tree = read_tree(path)
            matrix = rand(rng, 0:3, n, p)
            matrix[rand(rng, n, p) .< 0.7] .= 0
            matrix[end, :] .= 0
            matrix[:, end] .= 0
            samples = ["site_$i" for i in 1:n]
            row_order, col_order = randperm(rng, n), randperm(rng, p - 1)
            table = DataFrame(sample = String[], species = String[], value = Int[])
            # Split every abundance into duplicate records and shuffle records;
            # explicit domains determine the output order and retain empty rows.
            for i in 1:n, j in 1:p-1
                value = matrix[i, j]
                value == 0 && continue
                push!(table, (samples[i], tree.tip_labels[j], 1))
                push!(table, (samples[i], tree.tip_labels[j], value - 1))
            end
            table = table[randperm(rng, nrow(table)), :]
            data = CommunityData(table; sample_labels = samples[row_order],
                species_labels = tree.tip_labels[col_order])
            @test data.matrix == matrix[row_order, col_order]
            padded = matrix[row_order, :]
            @test align_community(data, tree).matrix == padded
            oracle = community_reference(padded, tree)
            @test same_values(pd(data, tree).pd, oracle.pd)
            @test same_values(phylo_endemism(data, tree).pe, oracle.pe)
            @test same_values(relative_pd(data, tree).rpd, oracle.rpd)
            @test same_values(beta_diss(data).total, taxonomic_beta_reference(padded, :sorensen).total)
            @test same_values(phylobeta(data, tree).total, oracle.beta(:sorensen).total)
            @test data.sample_labels == samples[row_order]
            @test pd(data, tree).sample_labels == samples[row_order]
        end
    end
end
