using Test, EvoKits, Random, SparseArrays, CSV, DataFrames

function timeslice_test_tree(newick)
    return mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, newick)
        read_tree(path)
    end
end

# Independent oracle: walk each ancestor path in BigFloat, then clip intervals
# in coordinates measured forward from the root. No descendant-age traversal.
function timeslice_reference(tree, depth)
    incoming = Dict(Int(tree.child_of_edge[e]) => BigFloat(tree.edge_length[e]) for e in 1:tree.nedges)
    depths = zeros(BigFloat, tree.nnodes)
    for node in 1:tree.nnodes
        current = node
        while current != tree.root
            depths[node] += incoming[current]
            current = Int(tree.parent_of_node[current])
        end
    end
    goal = max(0, maximum(depths[tree.tip_ids]) - BigFloat(depth))
    clipped = min.(depths, goal)
    lengths = [clipped[Int(tree.child_of_edge[e])] - clipped[Int(tree.parent_of_edge[e])]
        for e in 1:tree.nedges]
    return Float64.(lengths), Float64.(clipped)
end

function timeslice_reference_metrics(tree, comm, lengths)
    present = [any(comm[s, t] > 0 for t in descendant_tips(tree, tree.child_of_edge[e]))
        for s in axes(comm, 1), e in 1:tree.nedges]
    occupancy = vec(sum(present; dims = 1))
    pd_values = [sum(BigFloat(lengths[e]) for e in 1:tree.nedges if present[s, e]; init = BigFloat(0))
        for s in axes(comm, 1)]
    pe_values = [sum(BigFloat(lengths[e]) / occupancy[e] for e in 1:tree.nedges if present[s, e]; init = BigFloat(0))
        for s in axes(comm, 1)]
    strict = [sum(BigFloat(lengths[e]) for e in 1:tree.nedges if present[s, e] && occupancy[e] == 1;
        init = BigFloat(0)) for s in axes(comm, 1)]
    return Float64.(pd_values), Float64.(pe_values), Float64.(strict)
end

@testset "Time slices: hand calculations and exact cut boundaries" begin
    tree = timeslice_test_tree("(((A:1,B:1)AB:1,C:2)ABC:2,(D:2,E:2)DE:2)R;")
    # The Newick reader omits root labels; attach metadata explicitly here.
    tree.node_labels[Int(tree.root)] = "R"
    @test timeslice(tree, 0).edge_length == tree.edge_length
    @test timeslice(tree, 1.5).edge_length == [2, 2, 0.5, 0.5, 0, 0, 0.5, 0.5]
    @test timeslice(tree, 2).edge_length == [2, 2, 0, 0, 0, 0, 0, 0]
    @test timeslice(tree, 3).edge_length == [1, 1, 0, 0, 0, 0, 0, 0]
    @test all(iszero, timeslice(tree, 4).edge_length)
    @test all(==(2), timeslice(tree, 2).dist_from_root[tree.tip_ids])
    @test all(==(0), timeslice(tree, 4).dist_from_root)
    @test any(>(0), timeslice(tree, prevfloat(2.0)).edge_length[3:end])
    @test all(iszero, timeslice(tree, nextfloat(2.0)).edge_length[3:end])
    @test timeslice(timeslice(tree, 0.75), 1.25).edge_length == timeslice(tree, 2).edge_length
    @test timeslice(tree, 2; atol = 0.1).edge_length == timeslice(tree, 2).edge_length
    for depth in (0, 1, 1.5, 2, 4)
        sliced = timeslice(tree, depth)
        @test sliced.tip_ids == tree.tip_ids && sliced.tip_labels == tree.tip_labels
        @test sliced.node_labels == tree.node_labels
        @test sliced.parent_of_edge == tree.parent_of_edge && sliced.child_of_edge == tree.child_of_edge
        @test sliced.preorder == tree.preorder && sliced.postorder == tree.postorder
    end
    collapsed = timeslice(tree, 2; collapse = true)
    @test collapsed.nnodes == tree.nnodes - 1
    @test collapsed.tip_labels == tree.tip_labels
    @test "AB" ∉ collapsed.node_labels && "ABC" in collapsed.node_labels && "DE" in collapsed.node_labels
    @test collapsed.node_labels[Int(collapsed.root)] == "R"
    @test pairwise_distance(collapsed) == pairwise_distance(timeslice(tree, 2))
    star = timeslice(tree, 4; collapse = true)
    @test star.nnodes == tree.ntips + 1 && star.nedges == tree.ntips
    @test all(iszero, star.edge_length) && star.tip_labels == tree.tip_labels
    @test length(star.children[Int(star.root)]) == tree.ntips

    zero = timeslice_test_tree("((A:0,B:0):0,C:0);")
    @test timeslice(zero, 0; collapse = true).ntips == 3
    @test timeslice(zero, 0; collapse = true).nnodes == 4
    unary = timeslice_test_tree("((((A:1,B:1)AB:0)U:0)V:1,C:2)R;")
    @test timeslice(unary, 0; collapse = true).nnodes == 5
    @test timeslice(unary, 0; collapse = true).tip_labels == unary.tip_labels
    @test pairwise_distance(timeslice(unary, 0; collapse = true)) == pairwise_distance(unary)
    singleton = timeslice_test_tree("(A:2)R;")
    @test timeslice(singleton, 1; collapse = true).edge_length == [1]
    @test timeslice(singleton, 2; collapse = true).ntips == 1
    @test timeslice(timeslice_test_tree("(A:0);"), 0).ntips == 1
    tiny = timeslice_test_tree("((A:1,B:1)AB:1e-12,C:1.000000000001)R;")
    @test timeslice(tiny, 0; collapse = true).nnodes == tiny.nnodes
    for scale in (1e-200, 1e200)
        scaled = deepcopy(tree)
        scaled.edge_length .*= scale
        @test timeslice(scaled, 1.5scale).edge_length ≈ timeslice(tree, 1.5).edge_length .* scale
    end
    stem = timeslice_test_tree("(((A:1,B:1):1,C:2):1e20);")
    @test timeslice(stem, 1).edge_length == [1e20, 1, 1, 0, 0]
    # A stale distance cache does not affect either slicing or rebuilt caches.
    stale = deepcopy(tree)
    stale.dist_from_root .= NaN
    @test timeslice(stale, 2).edge_length == timeslice(tree, 2).edge_length
    @test timeslice(stale, 2).dist_from_root == timeslice(tree, 2).dist_from_root
    @test timeslice(stale, 2; collapse = true).dist_from_root == collapsed.dist_from_root
end

@testset "Time slices: independent path clipping and local R fixtures" begin
    path = joinpath(@__DIR__, "fixtures", "lineages")
    reference = CSV.read(joinpath(path, "timeslice_r.tsv"), DataFrame; delim = '\t')
    trees = unique(reference[:, [:case, :input_newick]])
    rng = MersenneTwister(20261107)
    for row in eachrow(trees)
        tree = timeslice_test_tree(row.input_newick)
        height = get_clades(tree; k = 1).tree_height
        cuts = unique(vcat(height .* [0, 0.13, 0.37, 0.71, 1],
            clamp.(height .- tree.dist_from_root[.!tree.is_tip], 0, height)))
        comm = rand(rng, 0:2, 5, tree.ntips)
        comm[end, :] .= 0
        for depth in cuts
            expected_lengths, expected_depths = timeslice_reference(tree, depth)
            plain = timeslice(tree, depth)
            collapsed = timeslice(tree, depth; collapse = true)
            @test plain.edge_length ≈ expected_lengths atol = 2e-14
            @test plain.dist_from_root ≈ expected_depths atol = 2e-14
            @test pairwise_distance(collapsed) ≈ pairwise_distance(plain) atol = 2e-14
            @test collapsed.tip_labels == tree.tip_labels
            @test all(collapsed.is_tip[collapsed.child_of_edge[e]] || collapsed.edge_length[e] > 0
                for e in 1:collapsed.nedges)
            p, pe, strict = timeslice_reference_metrics(tree, comm, expected_lengths)
            for sliced in (plain, collapsed)
                @test pd(comm, sliced).pd ≈ p atol = 2e-14
                @test phylo_endemism(comm, sliced).pe ≈ pe atol = 2e-14
                @test phylo_endemism(comm, sliced; weighted = false).pe ≈ strict atol = 2e-14
            end
        end
        for saved in eachrow(reference[reference.case .== row.case, :])
            actual = timeslice(tree, saved.depth; collapse = Bool(saved.collapse))
            expected = timeslice_test_tree(saved.newick)
            @test actual.tip_labels == expected.tip_labels
            @test actual.nnodes == expected.nnodes
            # Compare topology by descendant sets, independent of R/Julia IDs.
            actual_groups = sort([join(sort(descendant_tips(actual, node; labels = true)), "|")
                for node in actual.postorder_internal])
            expected_groups = sort([join(sort(descendant_tips(expected, node; labels = true)), "|")
                for node in expected.postorder_internal])
            @test actual_groups == expected_groups
            @test pairwise_distance(actual) ≈ pairwise_distance(expected) atol = 3e-14
            @test actual.dist_from_root[actual.tip_ids] ≈ expected.dist_from_root[expected.tip_ids] atol = 3e-14
            @test pd(comm, actual).pd ≈ pd(comm, expected).pd atol = 3e-14
        end
    end
end

@testset "Time slices: sparse community, beta and regionalization integration" begin
    tree = timeslice_test_tree("(((A:1,B:1):1,C:2):2,(D:2,E:2):2);")
    matrix = [1 0 1 0 0; 0 1 0 1 0; 0 0 0 0 0; 0 0 0 0 1]
    names = ["west", "mixed", "empty", "east"]
    data = CommunityData(sparse(matrix); species_labels = tree.tip_labels, sample_labels = names)
    ctx = prepare_community(data, tree)
    original_pd = copy(pd(ctx).pd)
    @test pd(data, timeslice(tree, 2)).pd == [2, 4, 0, 2]
    @test phylo_endemism(data, timeslice(tree, 2)).pe == [1, 2, 0, 1]
    for depth in (0, 0.5, 1, 1.5, 2, 3, 4)
        plain = timeslice(tree, depth)
        collapsed = timeslice(tree, depth; collapse = true)
        prepared = prepare_community(data, collapsed)
        @test pd(prepared).pd ≈ pd(data, plain).pd
        @test phylo_endemism(prepared).pe ≈ phylo_endemism(data, plain).pe
        @test pd(prepared).sample_labels == names && pd(prepared).richness == vec(sum(matrix; dims = 2))
        @test all(pd(prepared).pd .<= original_pd)
        for family in (:sorensen, :jaccard)
            beta = phylobeta(data, plain; index_family = family)
            other = phylobeta(prepared; index_family = family)
            for field in (:total, :turnover, :nestedness)
                @test isapprox(getfield(beta, field), getfield(other, field); nans = true)
            end
            selected = phylobeta(data, collapsed; index_family = family, pairs = [(1, 4), (2, 1), (3, 3)])
            @test isequal(selected.total, [beta.total[1, 4], beta.total[2, 1], beta.total[3, 3]])
        end
        if depth < 4
            subset = CommunityData(data.matrix[[1, 2, 4], :]; species_labels = data.species_labels,
                sample_labels = names[[1, 2, 4]])
            first = phyloregion(phylobeta(subset, plain); k = 2)
            second = phyloregion(phylobeta(subset, collapsed); k = 2)
            @test first.membership == second.membership
            @test first.region_distances ≈ second.region_distances
            @test first.sample_labels == names[[1, 2, 4]]
        else
            @test_throws ArgumentError phyloregion(phylobeta(data, plain); k = 2)
        end
    end
    @test pd(ctx).pd == original_pd
    # Partial, reordered species inputs retain their label-based interpretation.
    order = [5, 1, 4, 2]
    partial = CommunityData(sparse(matrix[:, order]); species_labels = tree.tip_labels[order], sample_labels = names)
    expected = copy(matrix)
    expected[:, 3] .= 0
    sliced = timeslice(tree, 1.5; collapse = true)
    @test phylo_endemism(partial, sliced).pe == phylo_endemism(expected, sliced).pe
    @test partial.matrix isa SparseMatrixCSC && partial.matrix == matrix[:, order]
    @test isequal(phylobeta(CommunityData(spzeros(Int, 0, 5); species_labels = tree.tip_labels), sliced).total,
        zeros(0, 0))
end

@testset "Time slices: validation and independent tree ownership" begin
    tree = timeslice_test_tree("((A:1,B:1):1,C:2);")
    for depth in (-1, Inf, NaN, true, "1", 3, BigFloat("1e-1000"))
        @test_throws ArgumentError timeslice(tree, depth)
    end
    for tolerance in (-1, Inf, NaN, true)
        @test_throws ArgumentError timeslice(tree, 1; atol = tolerance)
        @test_throws ArgumentError timeslice(tree, 1; rtol = tolerance)
    end
    @test_throws ArgumentError timeslice(timeslice_test_tree("(A:1,B:2);"), 0)
    near = timeslice_test_tree("(A:1,B:1.000001);")
    @test_throws ArgumentError timeslice(near, 0)
    @test timeslice(near, 0; atol = 1e-5).edge_length == near.edge_length
    for value in (-1.0, Inf, NaN)
        invalid = deepcopy(tree)
        invalid.edge_length[1] = value
        @test_throws ArgumentError timeslice(invalid, 0)
    end
    invalid = deepcopy(tree)
    invalid.edge_length .= 1e308
    @test_throws ArgumentError timeslice(invalid, 0)
    for label in ("", "  ", "B")
        invalid = deepcopy(tree)
        invalid.tip_labels[1] = label
        @test_throws ArgumentError timeslice(invalid, 0)
    end
    before = deepcopy(tree)
    for collapse in (false, true)
        sliced = timeslice(tree, 1; collapse)
        for field in fieldnames(CompactTree)
            original_value, new_value = getfield(tree, field), getfield(sliced, field)
            if original_value isa AbstractArray || original_value isa AbstractDict
                @test original_value !== new_value
            end
        end
        @test all(sliced.children[i] !== tree.children[i] for i in 1:min(sliced.nnodes, tree.nnodes))
        sliced.edge_length .= 99
        sliced.dist_from_root .= 99
        sliced.tip_labels[1] = "changed"
        empty!(sliced.children[Int(sliced.root)])
        @test all(isequal(getfield(tree, field), getfield(before, field)) for field in fieldnames(CompactTree))
    end
end
