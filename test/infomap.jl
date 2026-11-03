using Test, EvoKits, SparseArrays, Random, StableRNGs, LinearAlgebra, CSV, DataFrames
import IGraphs

# Independent dense two-level map equation for an undirected graph, without
# regularization. The test gets node assignments, not a backend-computed score.
function infomap_test_codelength(community, membership; weighted=false)
    columns = findall(vec(any(>(0), community; dims=1)))
    x = BigFloat.(community[:, columns])
    weighted || (x = BigFloat.(x .> 0))
    n, p = size(x)
    adjacency = [zeros(BigFloat, n, n) x; transpose(x) zeros(BigFloat, p, p)]
    flow = adjacency / sum(adjacency)
    stationary = vec(sum(flow; dims=2))
    plogp(x) = iszero(x) ? BigFloat(0) : x * log2(x)
    exits, usage = BigFloat[], BigFloat[]
    for region in unique(membership)
        inside = membership .== region
        exit = sum(flow[inside, .!inside])
        push!(exits, exit)
        push!(usage, exit + sum(stationary[inside]))
    end
    return Float64(plogp(sum(exits)) - 2sum(plogp, exits) + sum(plogp, usage) - sum(plogp, stationary))
end

@testset "Infomap: regularized defaults and independent map equation" begin
    blocks = [1 1 0 0; 1 1 0 0; 0 0 1 1; 0 0 1 1]
    plain = infomap(blocks; regularized=false, seed=42)
    regularized = infomap(blocks; seed=42)
    weak = infomap(blocks; regularization_strength=0.1, seed=42)
    @test plain.membership == [1, 1, 2, 2]
    @test plain.cluster_sizes == [2, 2] && plain.k == 2
    @test plain.codelength ≈ 2.0
    @test regularized.regularized && regularized.regularization_strength == 1.0
    @test regularized.membership == ones(Int, 4)
    @test regularized.cluster_sizes == [4] && regularized.k == 1
    @test regularized.codelength ≈ 3.0
    @test weak.membership == plain.membership
    @test weak.codelength ≈ plain.codelength
    @test regularized.seed == UInt64(42) && regularized.trials == 10
    @test regularized.weights === :binary && !isempty(regularized.igraph_version)
    @test regularized.sample_labels == ["sample_$i" for i in 1:4]
    for reg in (false, true)
        single = infomap(reshape([2.0], 1, 1); regularized=reg, weights=:abundance)
        @test single.membership == [1] && single.cluster_sizes == [1] && single.k == 1
        @test single.codelength ≈ 1.0
    end
    rng = StableRNG(907)
    for case in 1:8, mode in (:binary, :abundance)
        x = rand(rng, 0:4, 3 + case % 4, 3 + case % 3)
        x[:, 1] .+= 1
        network = EvoKits._infomap_network(x, mode)
        nodes, code, _ = EvoKits._infomap_backend(network.edges, network.edge_weights,
            network.nvertices, 5, UInt64(42), false, 1.0)
        reference = infomap_test_codelength(x, nodes; weighted=mode === :abundance)
        @test code ≈ reference atol=1e-12 rtol=1e-12
        fit = infomap(x; weights=mode, regularized=false, trials=5, seed=42)
        @test fit.codelength == code
        @test fit.membership == EvoKits._infomap_site_membership(nodes, size(x, 1))[1]
    end
end

@testset "Infomap: dense/CSC parity, labels, weights, zeros and ownership" begin
    x = [8.0 2 0 0; 1 9 0 0; 0 0 3 7; 0 0 9 1]
    samples = ["west_1", "west_2", "east_1", "east_2"]
    species = copy(samples) # Names may overlap across the two node types.
    before = copy(x)
    for mode in (:binary, :abundance), reg in (false, true), strength in (0.1, 1.0)
        a = infomap(x; sample_labels=samples, species_labels=species,
            weights=mode, regularized=reg, regularization_strength=strength, seed=89)
        b = infomap(sparse(x); sample_labels=samples, species_labels=species,
            weights=mode, regularized=reg, regularization_strength=strength, seed=89)
        data = CommunityData(x; sample_labels=samples, species_labels=species)
        c = infomap(data; weights=mode, regularized=reg, regularization_strength=strength, seed=89)
        @test a.membership == b.membership == c.membership
        @test a.codelength == b.codelength == c.codelength
        @test a.cluster_sizes == b.cluster_sizes == c.cluster_sizes
        @test a.sample_labels == samples && a.sample_labels !== samples
        @test all(1 .<= a.membership .<= a.k) && sum(a.cluster_sizes) == 4
        @test sort(unique(a.membership)) == collect(1:a.k)
        @test [findfirst(==(k), a.membership) for k in 1:a.k] |> issorted
        @test x == before && data.matrix == x
    end
    padded = sparse(hcat(zeros(4), x, zeros(4)))
    padded[1, 1] = 1.0
    padded[1, 1] = 0.0 # Retain an explicitly stored zero.
    @test nnz(padded) > count(>(0), padded)
    actual = infomap(padded; seed=42)
    @test actual.membership == infomap(x; seed=42).membership
    @test actual.codelength == infomap(x; seed=42).codelength
    @test nnz(padded) > count(>(0), padded) # Input is not cleaned in place.
    network = EvoKits._infomap_network(padded, :abundance)
    @test network.active_species == [2, 3, 4, 5]
    @test network.nvertices == 8
    @test network.edges == [0, 4, 1, 4, 0, 5, 1, 5, 2, 6, 3, 6, 2, 7, 3, 7]
    @test network.edge_weights == [8, 1, 2, 9, 3, 9, 7, 1]
    binary = infomap(x; regularized=false)
    @test binary.membership == infomap(x .> 0; regularized=false).membership
    @test binary.codelength == infomap(x .> 0; regularized=false).codelength
    huge = fill(big"1e1000", 2, 1)
    @test isfinite(infomap(huge).codelength) # Binary presence needs no Float64 conversion.
    @test_throws ArgumentError infomap(huge; weights=:abundance)
    @test_throws ArgumentError infomap(fill(big"1e-1000", 2, 1); weights=:abundance)
    @test_throws ArgumentError infomap(fill(1e308, 2, 1); weights=:abundance)
    fit = infomap(x; sample_labels=[1, 2, 3, 4])
    @test fit.sample_labels == ["1", "2", "3", "4"]
    fit.sample_labels[1] = "changed"
    fit.membership[1] = 99
    @test samples[1] == "west_1"
    @test infomap(x).membership != fit.membership
end

@testset "Infomap: reproducibility and independent native/Julia random streams" begin
    x = [1 1 0 0; 1 1 1 0; 0 1 1 1; 0 0 1 1; 1 0 0 1]
    for reg in (true, false), seed in (0, 42, typemax(UInt64))
        a = infomap(x; regularized=reg, seed, trials=5)
        b = infomap(x; regularized=reg, seed, trials=5)
        @test a.membership == b.membership && a.codelength == b.codelength
    end
    first = infomap(x; regularized=false, trials=1, seed=7)
    many = infomap(x; regularized=false, trials=10, seed=7)
    @test many.codelength <= first.codelength
    Random.seed!(603)
    expected = rand(UInt64)
    Random.seed!(603)
    infomap(x)
    @test rand(UInt64) == expected

    lib = IGraphs.LibIGraph
    generator = Ref{IGraphs.LibIGraph.igraph_rng_t}()
    control = Ref{IGraphs.LibIGraph.igraph_rng_t}()
    kind = cglobal((:igraph_rngtype_pcg32, IGraphs.LibIGraph.libigraph),
        IGraphs.LibIGraph.igraph_rng_type_t)
    @test lib.igraph_rng_init(generator, kind) == lib.IGRAPH_SUCCESS
    @test lib.igraph_rng_init(control, kind) == lib.IGRAPH_SUCCESS
    try
        @test lib.igraph_rng_seed(generator, UInt64(99)) == lib.IGRAPH_SUCCESS
        @test lib.igraph_rng_seed(control, UInt64(99)) == lib.IGRAPH_SUCCESS
        GC.@preserve generator control begin
            previous = lib.igraph_rng_set_default(generator)
            try
                infomap(x)
                @test lib.igraph_rng_default() == Base.unsafe_convert(Ptr{IGraphs.LibIGraph.igraph_rng_t}, generator)
                @test lib.igraph_rng_get_integer(lib.igraph_rng_default(), 0, 1000000) ==
                    lib.igraph_rng_get_integer(control, 0, 1000000)
                @test_throws ArgumentError infomap(x; trials=0)
                @test lib.igraph_rng_default() == Base.unsafe_convert(Ptr{IGraphs.LibIGraph.igraph_rng_t}, generator)
            finally
                lib.igraph_rng_set_default(previous)
            end
        end
    finally
        lib.igraph_rng_destroy(generator)
        lib.igraph_rng_destroy(control)
    end
end

@testset "Infomap: input validation, aggregation and sparse construction storage" begin
    x = [1 0; 0 1]
    for value in (zeros(0, 2), zeros(2, 0), zeros(2, 2), [1 0; 0 0],
                  [-1 0; 0 1], [NaN 0; 0 1], [Inf 0; 0 1], [missing 0; 0 1])
        @test_throws ArgumentError infomap(value)
        @test_throws ArgumentError infomap(sparse(value))
    end
    for trials in (0, -1, true, 1.5, big(typemax(Int)) + 1)
        @test_throws ArgumentError infomap(x; trials)
    end
    for seed in (-1, true, 0.5, big(typemax(UInt64)) + 1)
        @test_throws ArgumentError infomap(x; seed)
    end
    for strength in (0, -1, NaN, Inf, true, big"1e1000", big"1e-1000")
        @test_throws ArgumentError infomap(x; regularization_strength=strength)
    end
    @test_throws ArgumentError infomap(x; regularized=1)
    @test_throws ArgumentError infomap(x; weights=:count)
    @test_throws ArgumentError infomap(x; sample_labels=["same", "same"])
    @test_throws ArgumentError infomap(x; species_labels=["A"])
    @test_throws ArgumentError infomap(x; sample_labels=["", "B"])
    data = CommunityData(x)
    @test_throws ArgumentError infomap(data; species_labels=["A", "B"])
    @test_throws ArgumentError infomap(data; sample_labels=["a", "b"])
    @test_throws ArgumentError infomap(data; tip_labels=["A", "B"])
    grouped = aggregate_species(CommunityData([1 0 1 0; 0 1 0 1]), ["a", "a", "b", "b"])
    @test infomap(grouped).sample_labels == grouped.sample_labels
    @test isfinite(infomap(grouped).codelength)
    tree = mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
        read_tree(path)
    end
    community = CommunityData([1 1 0 0; 0 0 1 1]; species_labels=["A", "B", "C", "D"])
    lineages = collapse_range(community, tree; k=2)
    @test infomap(lineages.community; regularized=false).membership == [1, 2]

    n, p = 200, 100000
    large = sparse(repeat(1:n, 2), vcat(1:n, (p-n+1):p), ones(Int, 2n), n, p)
    EvoKits._infomap_network(large, :binary) # Warm compilation before allocation check.
    allocated = @allocated network = EvoKits._infomap_network(large, :binary)
    @test allocated < 4_000_000
    @test length(network.edges) == 4n && network.nvertices == 3n
    @test Base.summarysize(network) < 100000
end

@testset "Infomap: stable local phyloregion reference partitions" begin
    reference = CSV.read(joinpath(@__DIR__, "fixtures", "infomap", "reference.tsv"), DataFrame; delim='\t')
    for name in unique(reference.case)
        rows = reference[reference.case .== name, :]
        sort!(rows, :site)
        x = zeros(Int, first(rows.n), first(rows.p))
        for row in eachrow(rows)
            x[row.site, :] = parse.(Int, split(row.values, ','))
        end
        @test size(x, 1) == nrow(rows)
        for matrix in (x, sparse(x))
            fit = infomap(matrix; sample_labels=rows.label, regularized=false,
                seed=first(rows.seed), trials=first(rows.trials))
            @test fit.sample_labels == rows.label
            @test fit.membership == rows.membership
            @test fit.k == maximum(rows.membership)
            @test fit.cluster_sizes == [count(==(i), rows.membership) for i in 1:fit.k]
        end
    end
end
