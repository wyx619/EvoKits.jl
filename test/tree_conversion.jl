using Random, SparseArrays

# Recursive reference is deliberately confined to small trees. Its child
# order is the historical contract, including multifurcations and unary nodes.
function conversion_reference_order(tree)
    preorder, postorder = Int32[], Int32[]
    function visit(node)
        push!(preorder, node)
        for child in tree.children[node]
            visit(child)
        end
        push!(postorder, node)
    end
    visit(tree.root)
    return preorder, postorder
end

function conversion_test_newick(rng, p)
    parts = ["sp_$i:" * string(rand(rng, (0.0, 0.125, 1.0, 3.0))) for i in 1:p]
    while length(parts) > 1
        shuffle!(rng, parts)
        children = splice!(parts, 1:rand(rng, 2:min(5, length(parts))))
        parent = "(" * join(children, ",") * "):" * string(rand(rng, (0.0, 0.5, 2.0)))
        rand(rng, Bool) && (parent = "(" * parent * "):1")
        push!(parts, parent)
    end
    return only(parts) * ";"
end

function conversion_deep_newick(p, side)
    buffer = IOBuffer()
    if side === :left
        print(buffer, repeat("(", p - 1), "sp_1:1,sp_2:1)")
        for i in 3:p
            print(buffer, ":1,sp_", i, ":", i - 1, ")")
        end
    else
        for i in 1:p-2
            print(buffer, "(sp_", i, ":", p - i, ",")
        end
        print(buffer, "(sp_", p - 1, ":1,sp_", p, ":1)", repeat(":1)", p - 2))
    end
    print(buffer, ";")
    return String(take!(buffer))
end

@testset "Tree conversion preserves ordered traversal and geometry" begin
    rng = MersenneTwister(20261106)
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        examples = ["A;", "(A:0,B:1,C:2,D:3)root;", "((((A:1):2):3):4);",
            "(((A:1,B:2)ab:3,C:4)abc:5,(D:6,E:7,F:8)def:9)root;"]
        append!(examples, [conversion_test_newick(rng, rand(rng, 2:32)) for _ in 1:24])
        for newick in examples
            write_tree(path, newick)
            parsed = EvoKits.load_tree(path)
            tree = EvoKits.compact_tree(parsed)
            preorder, postorder = conversion_reference_order(tree)
            @test tree.preorder == preorder
            @test tree.postorder == postorder
            @test tree.postorder_internal == filter(node -> !tree.is_tip[node], postorder)
            @test sort(tree.postorder) == collect(1:tree.nnodes)
            @test tree.postorder[end] == tree.root
            @test all(tree.dist_from_root[tree.child_of_edge[e]] ==
                tree.dist_from_root[tree.parent_of_edge[e]] + tree.edge_length[e] for e in 1:tree.nedges)
            for other in (read_tree(path), EvoKits.preprocess_tree(parsed))
                @test all(isequal(getfield(other, field), getfield(tree, field)) for field in fieldnames(CompactTree))
            end
            if occursin("abc", newick)
                root_label = tree.node_labels[tree.root]
                @test tree.node_labels[tree.postorder] == ["A", "B", "ab", "C", "abc", "D", "E", "F", "def", root_label]
                @test tree.node_labels[tree.postorder_internal] == ["ab", "abc", "def", root_label]
            end
        end
    end
end

@testset "Tree conversion reads 65536-tip deep Newick combs" begin
    mktempdir() do dir
        path = joinpath(dir, "deep.nwk")
        p = 65_536
        for side in (:left, :right)
            write_tree(path, conversion_deep_newick(p, side))
            tree = read_tree(path)
            @test tree.ntips == p && tree.nnodes == 2p - 1 && tree.nedges == 2p - 2
            @test tree.tip_labels == ["sp_$i" for i in 1:p]
            @test tree.preorder == collect(1:tree.nnodes)
            @test all(==(p - 1), tree.dist_from_root[tree.tip_ids])
            if side === :left
                expected = Int32[p]
                for i in 2:p
                    push!(expected, p + i - 1, p - i + 1)
                end
                @test tree.postorder == expected
                @test tree.postorder_internal == collect(p-1:-1:1)
            else
                @test tree.postorder == [collect(2:2:2p-2); 2p - 1; collect(2p-3:-2:1)]
                @test tree.postorder_internal == collect(2p-3:-2:1)
            end
            @test distance(tree, "sp_1", "sp_$p") == 2(p - 1)
            # With patristic/2, the first and last tips are at distance p-1.
            # Two singleton sites therefore have zero alpha and pooled half that.
            input = sparse([1, 2], [1, p], [1.0, 1.0], 2, p)
            for pairwise in (false, true)
                result = raoD(input, tree; pairwise)
                @test result.Dkk == [0, 0]
                @test result.total == result.beta == (p - 1) / 2
                @test result.alpha == 0 && result.Fst == 1
                @test pairwise ? result.Dkl == [0 p-1; p-1 0] : size(result.Dkl) == size(result.H) == (0, 0)
            end
            if side === :left
                # The streaming Newick and NEXUS parsers share compact_tree.
                # Exercise both real input paths at the original failing depth.
                for stream_path in (path, joinpath(dir, "deep.nex"))
                    if stream_path != path
                        write(stream_path, "#NEXUS\nBEGIN TREES;\nTREE deep = " * read(path, String) * "\nEND;\n")
                    end
                    read_trees(stream_path) do stream
                        streamed = first(stream)
                        @test all(isequal(getfield(streamed, field), getfield(tree, field)) for field in fieldnames(CompactTree))
                        @test iterate(stream) === nothing
                    end
                end
            end
        end
        # Unary depth exercises the same conversion without relying on many tips.
        depth = 65_536
        write_tree(path, repeat("(", depth) * "only:1" * repeat("):1", depth) * ";")
        tree = read_tree(path)
        @test tree.ntips == 1 && tree.nnodes == depth + 1
        @test tree.tip_labels == ["only"] && only(tree.dist_from_root[tree.tip_ids]) == depth
        @test tree.postorder == collect(depth+1:-1:1)
        @test tree.postorder_internal == collect(depth:-1:1)
    end
end
