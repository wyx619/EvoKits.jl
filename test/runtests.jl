using Test
using EvoKits
using DataFrames
using Random

@testset "Tree file I/O" begin
    mktempdir() do dir
        input = joinpath(dir, "input.nwk")
        output = joinpath(dir, "output.nwk")
        write_tree(input, "((A:1,B:2)inner:3,C:4)root;")
        tree = read_tree(input)
        @test tree isa CompactTree
        loaded = EvoKits.compact_tree(EvoKits.load_tree(input))
        @test loaded.tip_labels == tree.tip_labels
        @test loaded.edge_length == tree.edge_length
        @test tree.tip_labels == ["A", "B", "C"]
        @test "inner" in tree.node_labels
        @test distance(tree, "A", "B") ≈ 3.0
        write_tree(output, tree)
        restored = read_tree(output)
        @test restored.tip_labels == tree.tip_labels
        @test restored.node_labels == tree.node_labels
        @test pairwise_distance(restored) ≈ pairwise_distance(tree)
        @test restored.edge_length == tree.edge_length
        write_tree(output, tree; append = true)
        @test length(readlines(output)) == 2
        @test_throws ArgumentError read_tree(joinpath(dir, "missing.nwk"))
        write(input, "(A:1,B:1)")
        @test_throws ArgumentError read_tree(input)
    end
end

@testset "EvoKits" begin
    tree = read_tree(joinpath(@__DIR__, "tiny.nwk"))
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
    tree = read_tree(tree_path)

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

@testset "EDGE2" begin
    tree_path = joinpath(mktempdir(), "edge2_tree.nwk")
    write(tree_path, "((A:1,B:1):1,(C:1,D:1):1);")
    tree = read_tree(tree_path)
    result = edge2_mod(tree, [0.1, 0.2, 0.3, 0.4])

    @test result.species == ["A", "B", "C", "D"]
    @test result.tbl == [1.0, 1.0, 1.0, 1.0]
    @test result.edge ≈ [0.12, 0.22, 0.42, 0.52]
    @test result.ed ≈ [1.2, 1.1, 1.4, 1.3]
    @test result.pd == 6.0
    @test result.epd_loss ≈ 1.14
    @test result.EDGE === result.edge
    @test result.PD == result.pd
    @test_throws ArgumentError edge2_mod(tree, [0.1, 0.2])
    @test_throws ArgumentError edge2_mod(tree, [0.1, 0.2, 0.0, 0.4])

    distribution = PextDistribution(
        ["LC", "NT", "VU", "EN", "CR"],
        [0.1, 0.2, 0.3, 0.4, 0.5],
    )
    sampled = sample_pext(["LC", "EX", "thr", "not", "unknown"], distribution; rng = MersenneTwister(2))
    @test sampled[2] == 1.0
    @test all((0.0 .< sampled) .& (sampled .<= 1.0))
    @test reorder_pext_to_tree(tree, ["D", "A", "B", "C"], [0.4, 0.1, 0.2, 0.3]) == [0.1, 0.2, 0.3, 0.4]
    @test_throws ArgumentError reorder_pext_to_tree(tree, ["A"], [0.1])
    adjusted = score_adjust([0.1, 0.3, 1.0], distribution, 0.2)
    @test length(adjusted) == 3
    @test adjusted[3] == 1.0
    @test all((0.0 .< adjusted) .& (adjusted .<= 1.0))

    csv_path = joinpath(mktempdir(), "pext.csv")
    write(csv_path, "RL.cat,pext\nLC,0.1\nNT,0.2\nVU,0.3\nEN,0.4\nCR,0.5\n")
    loaded = load_pext_distribution_csv(csv_path)
    @test loaded.categories == distribution.categories
    @test loaded.pext_values == distribution.pext_values

    generated = calc_pext_distribution(n_tips = 50; rng = MersenneTwister(3))
    @test length(generated.categories) == 50
    @test length(generated.pext_values) == 50
    @test all((0.0 .< generated.pext_values) .& (generated.pext_values .<= 1.0))
end

@testset "NRI and NTI local null model" begin
    tree = read_tree(joinpath(@__DIR__, "tiny.nwk"))
    mpd, mntd = EvoKits._community_metrics(tree, [1, 2, 3])
    @test mpd ≈ 10 / 3
    @test mntd ≈ 8 / 3

    result1 = nri_nti(tree, [1, 2, 3], ["A", "B"]; nreps = 30, seed = 41)
    result2 = nri_nti(tree, [1, 2, 3], [1, 2]; nreps = 30, seed = 41)
    @test result1.nri == result2.nri
    @test result1.nti == result2.nti
    @test result1.status == :ok
    @test result1.NRI == result1.nri
    @test result1.pool_richness == 3
    @test result1.observed_richness == 2
    @test stable_grid_seed(20260713, "grid-1") == stable_grid_seed(20260713, "grid-1")
    @test stable_grid_seed(20260713, "grid-1") != stable_grid_seed(20260713, "grid-2")

    degenerate = nri_nti(tree, [1, 2], [1, 2]; nreps = 10, seed = 1)
    @test degenerate.status == :null_sd_zero
    @test isnan(degenerate.nri)
    @test_throws ArgumentError nri_nti(tree, [1, 2], [1, 3]; nreps = 10)

    grid_seed = stable_grid_seed(20260713, "grid-1")
    grid_result = nri_nti(tree, [1, 2, 3], [1, 2]; nreps = 30,
                          seed = 20260713, grid_id = "grid-1")
    manual_seed_result = nri_nti(tree, [1, 2, 3], [1, 2]; nreps = 30,
                                  seed = grid_seed)
    @test grid_result.observed_mpd == manual_seed_result.observed_mpd
    @test grid_result.observed_mntd == manual_seed_result.observed_mntd
    @test grid_result.null_mpd_mean == manual_seed_result.null_mpd_mean
    @test grid_result.null_mntd_mean == manual_seed_result.null_mntd_mean
    @test grid_result.nri == manual_seed_result.nri
    @test grid_result.nti == manual_seed_result.nti
    @test grid_result.seed == grid_seed
    @test_throws ArgumentError nri_nti(tree, [1, 2, 3], [1, 2]; nreps = 30,
                                       rng = MersenneTwister(1), grid_id = "grid-1")

    table_dir = mktempdir()
    pool_table = joinpath(table_dir, "pool.tsv")
    observed_table = joinpath(table_dir, "observed.tsv")
    write(pool_table, "grid_id\tspecies\ngrid-a\tA\ngrid-a\tB\ngrid-a\tC\ngrid-b\tA\ngrid-b\tB\ngrid-b\tC\n")
    write(observed_table, "grid_id\tspecies\ngrid-a\tA\ngrid-a\tB\ngrid-b\tB\ngrid-b\tC\n")
    direct_csv = nri_nti(tree, ["A", "B", "C"], ["A", "B"];
                         nreps = 30, seed = 71, grid_id = "grid-a")
    table_serial = nri_nti_table(tree, pool_table, observed_table;
                                  nreps = 30, seed = 71, threaded = false, delim = '\t')
    table_threaded = nri_nti_table(tree, pool_table, observed_table;
                                    nreps = 30, seed = 71, threaded = true, delim = '\t')
    @test table_serial.grid_ids == ["grid-a", "grid-b"]
    @test table_serial.results[1].nri == direct_csv.nri
    @test table_serial.results[1].nti == direct_csv.nti
    @test table_threaded.results == table_serial.results

    auto_pool_table = joinpath(table_dir, "pool.txt")
    auto_observed_table = joinpath(table_dir, "observed.txt")
    write(auto_pool_table, "grid_id,species\ngrid-a,A\ngrid-a,B\ngrid-a,C\n")
    write(auto_observed_table, "grid_id,species\ngrid-a,A\ngrid-a,B\n")
    auto_table = nri_nti_table(tree, auto_pool_table, auto_observed_table;
                                nreps = 30, seed = 71, threaded = false)
    @test auto_table.results[1] == direct_csv

    noncontiguous_pool_table = joinpath(table_dir, "noncontiguous_pool.txt")
    write(noncontiguous_pool_table, "grid_id,species\ngrid-a,A\ngrid-b,A\ngrid-a,B\n")
    @test_throws ArgumentError nri_nti_table(tree, noncontiguous_pool_table, observed_table;
                                              nreps = 10, threaded = false, delim = ',')

    swq_like_path = joinpath(mktempdir(), "swq_like.nwk")
    write(swq_like_path,
          "(((A:0.1,B:0.2):1,(C:0.3,D:0.4):1):1,((E:0.5,F:0.6):1,(G:0.7,H:0.8):1):1);")
    swq_like_tree = read_tree(swq_like_path)
    query_index = EvoKits._build_tree_query_index(swq_like_tree)
    workspace = EvoKits._NriMetricsWorkspace(swq_like_tree)
    rng = MersenneTwister(20260729)
    for _ in 1:100
        selected = randperm(rng, swq_like_tree.ntips)[1:rand(rng, 2:swq_like_tree.ntips)]
        virtual_mpd, virtual_mntd = EvoKits._nri_community_metrics!(
            swq_like_tree, selected, query_index, workspace,
        )
        reference_mpd, reference_mntd = EvoKits._community_metrics(
            swq_like_tree, selected, query_index,
        )

        brute_mpd_total = 0.0
        brute_mntd_total = 0.0
        for i in eachindex(selected)
            nearest = Inf
            for j in eachindex(selected)
                i == j && continue
                value = distance(swq_like_tree, selected[i], selected[j])
                brute_mpd_total += i < j ? value : 0.0
                nearest = min(nearest, value)
            end
            brute_mntd_total += nearest
        end
        pair_count = length(selected) * (length(selected) - 1) / 2
        @test virtual_mpd ≈ brute_mpd_total / pair_count atol = 1e-12 rtol = 1e-12
        @test virtual_mntd ≈ brute_mntd_total / length(selected) atol = 1e-12 rtol = 1e-12
        @test virtual_mpd ≈ reference_mpd atol = 1e-12 rtol = 1e-12
        @test virtual_mntd ≈ reference_mntd atol = 1e-12 rtol = 1e-12
    end
end

@testset "Phylogenetic D" begin
    tree_path = joinpath(mktempdir(), "phylo_d_tree.nwk")
    write(tree_path, "((A:1,B:1):1,(C:1,D:1):1);")
    tree = read_tree(tree_path)

    result = phylo_d(tree, [1, 1, 0, 0]; permutations = 20, rng = MersenneTwister(11))
    @test result.observed ≈ 1.0
    @test result.d == result.D == result.DEstimate
    @test result.Pval1 == result.p_random
    @test result.Pval0 == result.p_brownian
    @test length(result.random_changes) == 20
    @test length(result.brownian_changes) == 20
    @test result.status == :ok

    stored = phylo_d(tree, [true, false, true, false]; permutations = 5, rng = MersenneTwister(12), store_nodal_values = true)
    @test size(stored.random_nodal_values) == (3, 5)
    @test size(stored.brownian_nodal_values) == (3, 5)
    @test stored.Parameters.Observed == stored.observed
    @test length(stored.NodalVals.observed) == 3

    df = DataFrame(species = ["D", "A", "C", "B"], state = [0, 1, 0, 1])
    aligned = phylo_d(tree, df, :state; taxon_col = :species, permutations = 5, rng = MersenneTwister(13))
    @test aligned.observed ≈ result.observed

    polytomy_path = joinpath(mktempdir(), "phylo_d_polytomy.nwk")
    write(polytomy_path, "(A:1,B:1,C:1);")
    polytomy = read_tree(polytomy_path)
    poly_result = phylo_d(polytomy, [1, 0, 0]; permutations = 5, rng = MersenneTwister(14))
    @test poly_result.observed ≈ 4 / 3

    subset = phylo_d_subset(tree, [1, 1, 0, 0]; permutations = 5, rng = MersenneTwister(15), min_tips = 2)
    @test length(subset.raw) == 3
    @test subset.tips == [2, 2, 4]
    @test subset.nodes == [1, 1, 3]
    @test subset.d[3] isa Float64

    @test_throws ArgumentError phylo_d(tree, [1, 1, 1, 1]; permutations = 5)
    @test_throws ArgumentError phylo_d(tree, [1, 1, 0]; permutations = 5)
    @test_throws ArgumentError phylo_d(tree, [1, 1, 0, 0]; permutations = 5, random_weights = [1.0, -1.0, 1.0, 1.0])
end

@testset "Tree kernel and shape" begin
    tree_path = joinpath(mktempdir(), "kernel_tree.nwk")
    write(tree_path, "((A:1,B:1):1,(C:1,D:1):1);")
    tree = read_tree(tree_path)

    @test mrca(tree, 1, 2) == Int(tree.parent_of_node[Int(tree.tip_ids[1])])
    @test mrca(tree, ["A", "B"]) == mrca(tree, 1, 2)
    @test distance(tree, 1, 2) == 2.0
    @test distance(tree, 1, 2; weighted = false) == 2
    @test pairwise_distance(tree) == [0.0 2.0 4.0 4.0; 2.0 0.0 4.0 4.0; 4.0 4.0 0.0 2.0; 4.0 4.0 2.0 0.0]
    @test distances_from_root(tree; weighted = false)[Int.(tree.tip_ids)] == [2, 2, 2, 2]
    @test distances_to_tip(tree, "A") == [0.0, 2.0, 4.0, 4.0]

    root = Int(tree.root)
    @test descendant_tips(tree, root; labels = true) == tree.tip_labels
    @test length(clade_members(tree)) == length(tree.postorder_internal)
    @test subtree_at_node(tree, root).tip_labels == tree.tip_labels

    @test tree_imbalance(tree; metric = :sackin) == 8.0
    @test tree_imbalance(tree; metric = :colless) == 0.0
    @test tree_imbalance(tree; metric = :colless_normalized) == 0.0
    @test tree_imbalance(tree; metric = :blum) ≈ log(3.0)
    @test gamma_statistic(tree) ≈ -0.8164965809277259
    ltt = lineages_through_time(tree; times = [0.0, 1.0, 2.0])
    @test ltt.lineages_at_query == [2, 4, 4]
end

@testset "P1 community beta and diversity" begin
    path = joinpath(mktempdir(), "star.nwk")
    write(path, "(A:1,B:1,C:1);")
    tree = read_tree(path)
    comm = [1 1 0; 0 1 1]

    @test comdist(comm, tree)[1, 2] == 1.5
    @test comdist(comm[:, [3, 1, 2]], tree; tip_labels = ["C", "A", "B"])[1, 2] == 1.5
    @test comdistnt(comm, tree)[1, 2] == 1.0
    @test unifrac(comm, tree)[1, 2] ≈ 2 / 3
    @test phylosor(comm, tree)[1, 2] ≈ 0.5
    rao = raoD(comm, tree)
    @test rao.Dkk == [0.5, 0.5]
    @test rao.Dkl[1, 2] == 0.75
    @test psv(comm, tree).psv ≈ [1.0, 1.0]
    @test psr(comm, tree).psr ≈ [2.0, 2.0]
    @test pse(comm, tree).pse ≈ [1.0, 1.0]
    @test psc(comm, tree).psc ≈ [1.0, 1.0]
    @test expected_pd(tree).expected_pd ≈ [1.0, 2.0, 3.0]
    @test edge_abundance_distribution(tree).num_children == [1, 1, 1]
    @test pd_calc(comm, tree; variant = :tip).pd == [2.0, 2.0]
    @test pd_calc(comm, tree; variant = :tbl).pd == [2.0, 2.0]
    @test pd_calc(comm, tree; variant = :sbl).pd == [0.0, 0.0]
    @test pd_calc(comm, tree; variant = :ueh).pd == [2.0, 2.0]
    boot = pd_bootstrap(comm, tree; reps = 3, rng = MersenneTwister(7))
    @test size(boot.replicates) == (2, 3)
    @test_throws ArgumentError comdist([1 -1 0], tree)

    balanced_path = joinpath(mktempdir(), "balanced_reference.nwk")
    write(balanced_path, "((A:1,B:1):1,(C:1,D:1):1);")
    balanced = read_tree(balanced_path)
    @test psv([1 1 0 0; 1 0 1 0; 1 1 1 1], balanced).variance ≈ [1 / 18, 1 / 18, 0.0]
    @test expected_pd(balanced).expected_pd ≈ [1.875, 3.5, 4.875, 6.0]
    @test expected_pd(balanced).variance_pd ≈ [3.2109375, 2.875, 1.3359375, 0.0]
    @test expected_pd(balanced; upper_bound = false).variance_pd ≈ [2.3671875, 2.375, 1.2421875, 0.0]
    @test edge_abundance_distribution(balanced).num_children == [2, 2, 1, 1, 1, 1]

    empty = [0 0 0 0]
    singleton = [1 0 0 0]
    @test pd(empty, balanced).pd == [0.0]
    @test pd(singleton, balanced).pd == [2.0]
    @test isnan(ses_mpd(singleton, balanced; runs = 3, rng = MersenneTwister(1)).mpd_obs[1])
    @test isnan(comdist(empty, balanced)[1, 1])
    @test unifrac(empty, balanced)[1, 1] == 0.0
    @test phylosor(empty, balanced)[1, 1] == 0.0
    @test isnan(psv(singleton, balanced).psv[1])
    nonultra_path = joinpath(mktempdir(), "nonultrametric.nwk")
    write(nonultra_path, "((A:1,B:2):1,(C:1,D:1):1);")
    nonultra = read_tree(nonultra_path)
    @test_throws ArgumentError raoD([1 0 0 0], nonultra)
end

@testset "Caper ED polytomy corrections" begin
    path = joinpath(mktempdir(), "ed_polytomy.nwk")
    write(path, "(A:1,B:1,C:1);")
    tree = read_tree(path)
    none = evol_distinct(tree; type = :fair_proportion, polytomy_correction = :none)
    isaac = evol_distinct(tree; type = :fair_proportion, polytomy_correction = :isaac)
    mooers = evol_distinct(tree; type = :fair_proportion, polytomy_correction = :mooers)
    @test none.scores == [1.0, 1.0, 1.0]
    @test isaac.scores ≈ fill(1.081 - 0.267 * log(3), 3)
    @test mooers.scores ≈ fill(0.8, 3)
    @test_throws ArgumentError evol_distinct(tree; polytomy_correction = :bad)
end

@testset "Picante equal-splits ordering" begin
    path = joinpath(mktempdir(), "ed_equal_splits.nwk")
    write(path, "(((A:1,B:1):3,C:1):2,D:1);")
    tree = read_tree(path)

    result = evol_distinct(tree; type = :equal_splits)
    @test result.tip_labels == ["A", "B", "C", "D"]
    @test result.scores ≈ [3.0, 3.0, 2.0, 1.0]
    @test evol_distinct(tree; type = :equal_splits, use_branch_lengths = false).scores ≈ [1.75, 1.75, 1.5, 1.0]
end

@testset "P2 restricted K and PIC" begin
    path = joinpath(mktempdir(), "pic.nwk")
    write(path, "((A:1,B:1):1,(C:1,D:1):1);")
    tree = read_tree(path)
    two_path = joinpath(mktempdir(), "two.nwk")
    write(two_path, "(A:1,B:1);")
    two = read_tree(two_path)
    pic_two = phylogenetic_contrasts(two, [1.0, 3.0])
    @test pic_two.contrasts ≈ [-sqrt(2.0)]
    @test pic_two.distances == [2.0]
    @test phylogenetic_contrasts(two, [1.0, 3.0]; scaled = false).contrasts == [-2.0]
    k = phylogenetic_signal_k(tree, [1.0, 2.0, 3.0, 4.0]; permutations = 20, rng = MersenneTwister(5))
    @test k.valid_replicates == 20
    @test isfinite(k.k)
    @test_throws ArgumentError phylogenetic_contrasts(tree, [1.0, 2.0, 3.0, NaN])
end
