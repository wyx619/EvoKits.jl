using Random, LinearAlgebra

function original_path_distance(tree, a, b)
    edges = Dict(tree.child_of_edge[e] => e for e in 1:tree.nedges)
    function path(node)
        result = Set{Int}()
        while node != tree.root
            edge = edges[node]
            push!(result, edge)
            node = tree.parent_of_node[node]
        end
        result
    end
    sum((BigFloat(tree.edge_length[e]) for e in symdiff(path(a), path(b))); init = BigFloat(0))
end

@testset "Public tree distances preserve original branch paths" begin
    mktempdir() do dir
        path = joinpath(dir, "distances.nwk")
        for newick in ("((A:1e-100,B:2e-100):1e200,C:1e200);",
            "((A:1000000001,B:1000000001):1e16,C:1);",
            "((A:1,B:1):1e308,C:1);", "((A:0,B:0):0,C:0);",
            "(((A:1,B:2):3,C:4):5,D:6);", "A;")
            write_tree(path, newick)
            tree = read_tree(path)
            all_nodes = collect(1:tree.nnodes)
            expected = [Float64(original_path_distance(tree, Int32(a), Int32(b))) for a in all_nodes, b in all_nodes]
            actual = pairwise_distance(tree; tips = all_nodes, by = :node)
            @test actual ≈ expected rtol=2e-14
            @test actual == transpose(actual)
            for a in all_nodes, b in all_nodes
                @test distance(tree, a, b; by = :node) ≈ expected[a, b] rtol=2e-14
            end
            for (i, tip) in enumerate(tree.tip_ids)
                @test distances_to_tip(tree, tree.tip_labels[i]) ≈ expected[tip, tree.tip_ids] rtol=2e-14
            end
            weighted = pairwise_distance(tree)
            @test weighted ≈ expected[tree.tip_ids, tree.tip_ids] rtol=2e-14
            @test diag(weighted) == zeros(tree.ntips)
            @test size(pairwise_distance(tree; tips = [])) == (0, 0)
            @test all(isinteger, pairwise_distance(tree; weighted = false))
            if occursin("1e-100", newick)
                @test distance(tree, "A", "B") ≈ 3e-100
                @test weighted[1, 2] ≈ 3e-100
            elseif occursin("1e308", newick)
                @test distance(tree, "A", "B") == 2
            end
        end
    end
end
