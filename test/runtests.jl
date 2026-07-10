using Test
using EvoKits
using DataFrames
using Random

@testset "EvoKits" begin
    nw = EvoKits.load_newick_tree(joinpath(@__DIR__, "tiny.nwk"))
    tree = EvoKits.to_compact_tree(nw)
    @test tree.ntips == 3

    ed = evol_distinct(tree)
    @test length(ed.scores) == 3

    ultra = make_ultrametric(tree)
    @test ultra.ntips == tree.ntips

    bif = multi2di(tree)
    @test bif.ntips == tree.ntips

    trait = [1.0, 2.0, 3.0]
    lam = pagel_lambda(ultra, trait)
    @test lam.success

    df = DataFrame(taxon = tree.tip_labels, x = trait, y = [1.0, missing, 3.0])
    aligned1 = align_traits_to_tree(tree, df; taxon_col = :taxon, trait_cols = :x)
    aligned2 = align_traits_to_tree(tree, df; taxon_col = :taxon, trait_cols = [:x, :y])
    @test length(aligned1) == 3
    @test size(aligned2) == (3, 2)

    comm = [1 0 1; 0 1 1]
    pdres = pd(comm, tree)
    @test length(pdres.pd) == 2

    mpdres = ses_mpd(comm, tree; runs = 10, rng = MersenneTwister(1))
    mntdres = ses_mntd(comm, tree; runs = 10, rng = MersenneTwister(1))
    spdres = ses_pd(comm, tree; runs = 10, rng = MersenneTwister(1))
    @test length(mpdres.mpd_obs) == 2
    @test length(mntdres.mntd_obs) == 2
    @test length(spdres.pd_obs) == 2
end

@testset "PhyloRef and tree pruning" begin
    tree_path = joinpath(mktempdir(), "mapping_tree.nwk")
    write(tree_path, "((A:1,B:1):1,(C:1,D:1):1);")
    tree = to_compact_tree(load_newick_tree(tree_path))

    map = build_phyloref(tree)
    @test map.ape_node_ids[Int(tree.root)] == tree.ntips + 1
    @test phyloref_node_table(tree; map = map).compact_node_id[Int(tree.root)] == Int(tree.root)
    edge_table = phyloref_edge_table(tree; map = map)
    @test "A|B" in edge_table.descendant_signature
    @test EvoKits.phylo_edges_from_signatures(tree, edge_table.descendant_signature; map = map) == collect(1:tree.nedges)
    @test ape_node_id(tree, Int(tree.root); map = map) == tree.ntips + 1
    @test compact_node_id(tree, tree.ntips + 1; map = map) == Int(tree.root)
    @test ape_edge_id(tree, 2; map = map) == 4
    @test compact_edge_id(tree, 4; map = map) == 2
    @test !(:preprocess_tree in names(EvoKits))
    @test ape_edge_matrix(tree; map = map) == [
        5 6;
        6 1;
        6 2;
        5 7;
        7 3;
        7 4
    ]
    node_attached = attach_ape_node_map(DataFrame(node_id = [Int(tree.root)]), tree; map = map)
    @test node_attached.ape_node_id == [tree.ntips + 1]
    edge_attached = attach_ape_edge_map(DataFrame(edge_id = [1]), tree; map = map)
    @test edge_attached.ape_edge_id == [5]

    pruned = keep_tip(tree, ["A", "C", "D"])
    @test Set(pruned.tip_labels) == Set(["A", "C", "D"])
    @test pruned.ntips == 3
    @test all(isapprox.(pruned.dist_from_root[pruned.tip_ids], 2.0))
    a_node = Int(pruned.tipname_to_id["A"])
    a_edge = findfirst(==(Int32(a_node)), pruned.child_of_edge)
    @test pruned.edge_length[a_edge] == 2.0

    dropped = drop_tip(tree, "B")
    @test Set(dropped.tip_labels) == Set(["A", "C", "D"])
    @test all(isapprox.(dropped.dist_from_root[dropped.tip_ids], 2.0))
    @test_throws ArgumentError keep_tip(tree, ["A", "missing"])
    @test_throws ArgumentError drop_tip(tree, tree.tip_labels)
end
