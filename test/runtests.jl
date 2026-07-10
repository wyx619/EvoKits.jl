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
