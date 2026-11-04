using SparseArrays

@testset "Community long-table records, ordering and ownership" begin
    matrix = Float32[0 2.5 0 0; 0 0 0 0; 3 0 0 0; 1e-5 4 0 0]
    labels = (; sample_labels = ["z", "empty", "群落甲", "a"],
        species_labels = ["C", "A", "absent", "B"])
    expected = DataFrame(sample = ["群落甲", "a", "z", "a"],
        species = ["C", "C", "A", "A"], value = Float32[3, 1e-5, 2.5, 4])
    for input in (matrix, view(matrix, :, :), sparse(matrix),
        SparseMatrixCSC{Float32, Int32}(sparse(matrix)))
        data = CommunityData(input; labels...)
        records = community_longtable(data)
        @test keys(records) == (:table, :sample_labels, :species_labels)
        @test isequal(records.table, expected)
        @test records.sample_labels == labels.sample_labels
        @test records.species_labels == labels.species_labels
        @test eltype(records.table.value) === Float32
        restored = CommunityData(records.table;
            sample_labels = records.sample_labels, species_labels = records.species_labels)
        @test restored.matrix isa SparseMatrixCSC{Float32}
        @test restored.matrix == matrix
        @test restored.sample_labels == labels.sample_labels
        @test restored.species_labels == labels.species_labels
        @test data.matrix == input
        # All outputs own their column/label arrays, and reconstruction copies them.
        records.table.value[1] = 99
        records.table.sample[1] = "changed"
        records.table.species[1] = "changed"
        records.sample_labels[1] = "changed"
        records.species_labels[1] = "changed"
        @test restored.matrix == data.matrix == matrix
        @test restored.sample_labels == data.sample_labels == labels.sample_labels
        @test restored.species_labels == data.species_labels == labels.species_labels
        fresh = community_longtable(data)
        data.matrix[1, 1] = 100
        data.sample_labels[1] = "mutated"
        data.species_labels[1] = "mutated"
        @test isequal(fresh.table, expected)
        @test fresh.sample_labels == labels.sample_labels
        @test fresh.species_labels == labels.species_labels
    end
    # Explicit zeros can appear after construction, and must not be exported
    # or removed from the caller's storage by a read-only conversion.
    data = CommunityData(sparse([1.0 0; 0 2]))
    nonzeros(data.matrix)[1] = 0
    records = community_longtable(data)
    @test records.table.value == [2.0]
    @test nnz(data.matrix) == 2
    @test nonzeros(data.matrix) == [0, 2]
end

@testset "Community long-table numeric values and empty label domains" begin
    for value in (true, Int8(3), typemax(Int), typemax(UInt), big(10)^100,
        Float32(0.125), nextfloat(0.0), 3 // 7, BigFloat("0.125"))
        T = typeof(value)
        for input in (reshape(T[value, 0, 0, 0], 2, 2),
            sparse([1], [1], T[value], 2, 2))
            data = CommunityData(input; sample_labels = [20, 10], species_labels = ["B", "A"])
            records = community_longtable(data)
            @test records.table.value == [value]
            @test eltype(records.table.value) === T
            restored = CommunityData(records.table;
                sample_labels = records.sample_labels, species_labels = records.species_labels)
            @test restored.matrix == input
            @test eltype(restored.matrix) === (T === Bool ? Int : T)
        end
        for shape in ((0, 0), (0, 2), (3, 0), (3, 2)), storage in (:preserve, :sparse)
            data = CommunityData(zeros(T, shape...); storage)
            records = community_longtable(data)
            @test nrow(records.table) == 0
            @test names(records.table) == ["sample", "species", "value"]
            @test eltype(records.table.value) === T
            restored = CommunityData(records.table;
                sample_labels = records.sample_labels, species_labels = records.species_labels)
            @test restored.matrix == data.matrix
            @test eltype(restored.matrix) === (T === Bool ? Int : T)
            @test restored.sample_labels == data.sample_labels
            @test restored.species_labels == data.species_labels
        end
    end
    # Empty untyped columns keep the constructor's Float64 fallback.
    empty = DataFrame(sample = String[], species = String[], value = Any[])
    @test eltype(CommunityData(empty).matrix) === Float64
    @test eltype(CommunityData(Int[], Int[], Real[];
        sample_labels = String[], species_labels = String[]).matrix) === Float64
end

@testset "Community long-table validation of modified snapshots" begin
    for storage in (:preserve, :sparse), bad in (-1.0, NaN, Inf)
        data = CommunityData(ones(2, 2); storage)
        data.matrix[1, 1] = bad
        @test_throws ArgumentError community_longtable(data)
    end
    for field in (:sample_labels, :species_labels)
        for invalid in (["a", "a"], ["", "b"], ["a"])
            data = CommunityData(ones(2, 2))
            labels = getproperty(data, field)
            empty!(labels)
            append!(labels, invalid)
            @test_throws ArgumentError community_longtable(data)
        end
    end
end

@testset "Community long-table randomized round trips and aggregation" begin
    rng = MersenneTwister(20261104)
    for trial in 1:25
        n, p = rand(rng, 2:12), rand(rng, 2:15)
        matrix = rand(rng, 0:8, n, p) ./ 4
        matrix[rand(rng, n, p) .< 0.7] .= 0
        matrix[end, :] .= 0
        matrix[:, end] .= 0
        labels = (; sample_labels = ["site_$i" for i in randperm(rng, n)],
            species_labels = ["sp_$j" for j in randperm(rng, p)])
        for storage in (:preserve, :sparse)
            data = CommunityData(matrix; storage, labels...)
            records = community_longtable(data)
            # Record order may change in external software; full labels restore
            # the exact domain and order without relying on first occurrences.
            shuffled = records.table[randperm(rng, nrow(records.table)), :]
            restored = CommunityData(shuffled; labels...)
            @test restored.matrix == matrix
            @test nrow(records.table) == count(>(0), matrix)
            @test isequal(beta_diss(restored).total, beta_diss(data).total)
        end
    end
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
        tree = read_tree(path)
        data = CommunityData(sparse([2.0 1 0 0; 0 0 3 1; 1 0 2 0; 0 0 0 0]);
            sample_labels = ["west", "east", "mixed", "empty"], species_labels = tree.tip_labels)
        groups = Dict("A" => "x", "B" => "x", "C" => "y", "D" => "y")
        traits = DataFrame(species = tree.tip_labels, size = [0.0, 0.0, 10.0, 10.0])
        grouping = functional_groups(data, traits; k = 2, nstart = 2, rng = MersenneTwister(10))
        outputs = (data, aggregate_species(data, groups; abundance_weighted = true, normalize = true),
            collapse_range(data, tree; cut = 1, statistic = :abundance).community,
            collapse_traits(data, grouping; statistic = :abundance, normalize = true))
        for output in outputs
            records = community_longtable(output)
            restored = CommunityData(records.table;
                sample_labels = records.sample_labels, species_labels = records.species_labels)
            @test restored.matrix == output.matrix
            @test restored.sample_labels == output.sample_labels
            @test restored.species_labels == output.species_labels
        end
        records = community_longtable(data)
        restored = CommunityData(records.table;
            sample_labels = records.sample_labels, species_labels = records.species_labels)
        @test pd(restored, tree).pd == pd(data, tree).pd
        @test phylo_endemism(restored, tree).pe == phylo_endemism(data, tree).pe
        @test isequal(phylobeta(restored, tree).total, phylobeta(data, tree).total)
    end
end

@testset "Community long-table sparse storage scales with records" begin
    n, p = 10_000, 50_000
    rows = [1, n, 42, 100]
    cols = [1, 2, p, p]
    data = CommunityData(rows, cols, [1.5, 2.5, 3.5, 1e100];
        sample_labels = ["site_$i" for i in 1:n], species_labels = ["sp_$j" for j in 1:p])
    exported = community_longtable(data)
    @test nrow(exported.table) == 4
    @test exported.table.value == [1.5, 2.5, 3.5, 1e100]
    # Dense storage would require 4 GB, and occurrence expansion is impossible.
    @test @allocated(community_longtable(data)) < 40_000_000
    @test Base.summarysize(exported) < 6_000_000
end
