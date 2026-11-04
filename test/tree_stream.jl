using Test, EvoKits, Random, StableRNGs, SparseArrays

function stream_test_files(f)
    return mktempdir() do dir
        first_path, second_path = joinpath(dir,"run1.nwk"), joinpath(dir,"run2.nwk")
        first_text = "\ufeff[&R] ((A:1,B:2):3,C:4);\n\n((C:3,A:2):1,B:5); ((A:2,C:4):1,B:3);\n"
        second_text = "[header]\n((B:3,C:2):2,A:1);\n((C:1,B:4):3,A:2);\n"
        write(first_path,first_text)
        write(second_path,second_text)
        return f(dir,[first_path,second_path])
    end
end

function same_stream_tree(left,right)
    return left.tip_labels == right.tip_labels && left.parent_of_edge == right.parent_of_edge &&
        left.child_of_edge == right.child_of_edge && isequal(left.edge_length,right.edge_length)
end

@testset "Tree streaming: records, files, burn-in, thinning and cleanup" begin
    stream_test_files() do dir,paths
        stream = read_trees(paths)
        @test stream isa TreeStream
        @test stream.io === nothing && isopen(stream)
        @test eltype(stream) === CompactTree
        trees = collect(stream)
        @test length(trees) == 5
        @test !isopen(stream) && stream.io === nothing
        @test iterate(stream) === nothing
        @test trees[1].tip_labels == ["A","B","C"]
        @test distance(trees[1],"A","B") == 3
        @test distance(trees[2],"A","C") == 5
        @test distance(trees[4],"B","C") == 5
        @test trees[2].tip_labels == ["C","A","B"]
        @test length(collect(read_trees(paths[1]))) == 3
        @test length(collect(read_trees((paths[2],paths[1])))) == 5
        kept = read_trees(paths;burnin=1,thin=2) do source
            collect(source)
        end
        @test length(kept) == 2
        @test same_stream_tree(kept[1],trees[2]) && same_stream_tree(kept[2],trees[5])
        @test isempty(collect(read_trees(paths;burnin=10)))
        empty_path = joinpath(dir,"empty.nwk")
        write(empty_path,"\n [only [nested] comments]\n")
        @test isempty(collect(read_trees(empty_path)))
        second_empty = empty_path * "2"
        write(second_empty,"")
        @test length(collect(read_trees([empty_path,paths[1],second_empty]))) == 3
    end
end

@testset "Tree streaming: comments, quoted labels and multiline Newick" begin
    mktempdir() do dir
        path = joinpath(dir,"quoted.nwk")
        write(path,"[&R] (('A;one':1[semi;colon],'B''quote':2)'inner,node':3,\n" *
            "'物种 C':4);\n('x[y]':2,'space name':3);\n")
        trees = collect(read_trees(path))
        @test length(trees) == 2
        @test trees[1].tip_labels == ["A;one","B'quote","物种 C"]
        @test "inner,node" in trees[1].node_labels
        @test distance(trees[1],"A;one","物种 C") == 8
        @test distance(trees[1],"A;one","B'quote") == 3
        @test trees[2].tip_labels == ["x[y]","space name"]
        @test distance(trees[2],"x[y]","space name") == 5
        @test trees[1].tipname_to_id["物种 C"] in trees[1].tip_ids
        write(path,"((A:1e-2,B:0):3,C:-2);\n((A,B),C);")
        trees = collect(read_trees(path))
        @test distance(trees[1],"A","B") ≈ 0.01
        @test any(==(-2),trees[1].edge_length)
        @test all(isnan,trees[2].edge_length)
        write(path,"((__EVOKITS_QUOTED_1:1,'quoted:tip':2):3,C:4);")
        @test collect(read_trees(path))[1].tip_labels == ["__EVOKITS_QUOTED_1","quoted:tip","C"]
        write(path,"single:2; 'single;quoted':1;")
        standalone = collect(read_trees(path))
        @test standalone[1].ntips == 1 && standalone[1].tip_labels == ["single"]
        @test standalone[2].tip_labels == ["single;quoted"]
    end
end

@testset "Tree streaming: lazy errors and deterministic closure" begin
    mktempdir() do dir
        path = joinpath(dir,"bad.nwk")
        write(path,"(A:1,B:1);\n(A:1,,B:1);")
        source = read_trees(path)
        @test first(source).ntips == 2
        handle = source.io
        observed_error = try
            iterate(source)
            nothing
        catch err
            err
        end
        @test observed_error isa ArgumentError
        @test occursin("tree 2, line 2",sprint(showerror,observed_error))
        @test occursin(abspath(path),sprint(showerror,observed_error))
        @test !isopen(source) && !isopen(handle)
        # Do-block closes the live stream even for an early return or callback error.
        write(path,"(A:1,B:1);(A:2,B:2);")
        captured = Ref{TreeStream}()
        tree = read_trees(path) do stream
            captured[] = stream
            first(stream)
        end
        @test tree.ntips == 2 && !isopen(captured[])
        @test_throws ErrorException read_trees(path) do stream
            captured[] = stream
            first(stream)
            error("callback stopped")
        end
        @test !isopen(captured[])
        stream = read_trees(path)
        first(stream)
        handle = stream.io
        close(stream)
        close(stream)
        @test !isopen(stream) && !isopen(handle)
        @test iterate(stream) === nothing
        # Syntax parsing is lazy; discarded records are framed but not parsed.
        write(path,"(A:1,,B:1);\n(A:1,B:2);")
        @test only(collect(read_trees(path;burnin=1))).ntips == 2
        for text in ("(A:1,B:1)","(A:1,'broken:1);","(A:1,B:1);[broken",
            "(A:1,B:1)]", "#NEXUS\nBEGIN TREES;", "(A:1,B:1));", ";", "();",
            "(A:1,);", "((A:1,B:1),);", "(A:bad,B:1);", "(A:1:2,B:1);",
            "(A:1,B:1) trailing text (C:1,D:1);", "(A:1,A:2);", "('':1,B:1);",
            "(A:1,'B':1'broken');")
            write(path,text)
            stream = read_trees(path)
            @test_throws ArgumentError collect(stream)
            @test !isopen(stream)
        end
        write(path,"(A:1,B:1)")
        @test_throws ArgumentError collect(read_trees(path;burnin=20))
        # Burn-in/thinning cannot turn a NEXUS container into Newick records.
        write(path,"#NEXUS\nBEGIN TREES;\nTREE one = (A:1,B:1);\nEND;")
        @test_throws ArgumentError collect(read_trees(path;burnin=2))
        @test_throws ArgumentError tree_sampler(path;n=1,burnin=2)
    end
end

@testset "Tree streaming: path and option validation" begin
    stream_test_files() do dir,paths
        for bad in ([],(),1,missing,[paths[1],1],(paths[1],paths[1]),
            [paths[1],joinpath(dir,".","run1.nwk")],joinpath(dir,"missing.nwk"),dir)
            @test_throws ArgumentError read_trees(bad)
        end
        for bad in (-1,true,1.,"1",missing,big(typemax(Int))+1)
            @test_throws ArgumentError read_trees(paths;burnin=bad)
            @test_throws ArgumentError read_trees(paths;thin=bad)
        end
        @test_throws ArgumentError read_trees(paths;thin=0)
    end
end

@testset "Posterior sampling: reproducibility, population weighting and provenance" begin
    stream_test_files() do dir,paths
        all_trees = collect(read_trees(paths))
        first = tree_sampler(paths;n=80,rng=StableRNG(42))
        second = tree_sampler(Tuple(paths);n=80,rng=StableRNG(42))
        @test first isa TreeSampleResult && eltype(first) === CompactTree
        @test length(first) == 80 && first.population_size == 5 && first.replace
        @test first.sample_indices == second.sample_indices
        expected_rng = StableRNG(42)
        @test first.sample_indices == [rand(expected_rng,1:5) for _ in 1:80]
        @test Set(first.sample_indices) == Set(1:5)
        @test all(same_stream_tree.(first.trees,all_trees[first.sample_indices]))
        @test all(same_stream_tree.(collect(first),collect(second)))
        for draw in 1:length(first)
            index = first.sample_indices[draw]
            @test first.source_paths[draw] == abspath(paths[index <= 3 ? 1 : 2])
            @test first.source_indices[draw] == (index <= 3 ? index : index-3)
            @test first.source_lines[draw] == [1,3,3,2,3][index]
        end
        duplicate = findfirst(i -> first.sample_indices[i] == first.sample_indices[1],2:length(first)) + 1
        saved = copy(first.trees[duplicate].edge_length)
        fill!(first.trees[1].edge_length,-1)
        @test first.trees[duplicate].edge_length == saved
        @test first.trees[1].tip_labels !== first.trees[duplicate].tip_labels
        without = tree_sampler(paths;n=5,replace=false,rng=StableRNG(7))
        @test !without.replace && length(unique(without.sample_indices)) == 5
        @test sort(without.sample_indices) == 1:5
        @test all(same_stream_tree.(without.trees,all_trees[without.sample_indices]))
        @test tree_sampler(paths;n=3,replace=false,rng=StableRNG(8)).sample_indices ==
            tree_sampler(paths;n=3,replace=false,rng=StableRNG(8)).sample_indices
        kept = tree_sampler(paths;n=2,replace=false,burnin=1,thin=2,rng=StableRNG(5))
        @test kept.population_size == 2 && kept.burnin == 1 && kept.thin == 2
        @test all(==(2),kept.source_indices)
        @test sort(kept.source_paths) == sort(abspath.(paths))
        @test all(same_stream_tree.(kept.trees,[all_trees[[2,5][i]] for i in kept.sample_indices]))
    end
end

@testset "Posterior sampling: invalid counts, populations and selected-only parsing" begin
    stream_test_files() do dir,paths
        for bad in (0,-1,true,1.,"1",missing,big(typemax(Int))+1)
            @test_throws ArgumentError tree_sampler(paths;n=bad)
        end
        @test_throws ArgumentError tree_sampler(paths;n=6,replace=false)
        @test_throws ArgumentError tree_sampler(paths;n=1,burnin=20)
        empty_path = joinpath(dir,"empty.nwk")
        write(empty_path," [comment] \n")
        @test_throws ArgumentError tree_sampler(empty_path;n=1)
        @test tree_sampler(paths[1];n=10).population_size == 3
        path = joinpath(dir,"sample.nwk")
        write(path,"(A:1,,B:1);\n(A:1,B:2);")
        @test only(tree_sampler(path;n=1,burnin=1)).ntips == 2
        # Find a deterministic seed selecting the syntactically valid record.
        seed = findfirst(seed -> rand(StableRNG(seed),1:2) == 2,1:20)
        sample = tree_sampler(path;n=1,rng=StableRNG(seed))
        @test sample.source_indices == [2]
        @test_throws ArgumentError tree_sampler(path;n=2,replace=false)
        write(path,"(A:1,B:1);\n(A:1,B:2)")
        @test_throws ArgumentError tree_sampler(path;n=1)
    end
end

@testset "Posterior sampling: distinct-draw distribution and bounded metadata" begin
    rng = StableRNG(20261104)
    frequencies = zeros(Int,5,5)
    for _ in 1:3000
        chosen = EvoKits._tree_sample_indices(rng,5,2,false)
        frequencies[chosen[1],chosen[2]] += 1
    end
    @test all(iszero(frequencies[i,i]) for i in 1:5)
    @test all(100 <= frequencies[i,j] <= 210 for i in 1:5 for j in 1:5 if i != j)
    EvoKits._tree_sample_indices(rng,10^12,8,false)
    bytes = @allocated chosen = EvoKits._tree_sample_indices(rng,10^12,8,false)
    @test length(unique(chosen)) == 8 && all(i -> 1 <= i <= 10^12,chosen)
    @test bytes < 20_000
    mktempdir() do dir
        path = joinpath(dir,"many.nwk")
        open(path,"w") do io
            for _ in 1:5000
                println(io,"((A:1,B:2):3,C:4);")
            end
        end
        stream = read_trees(path)
        try
            for _ in 1:2000
                iterate(stream)
            end
            @test Base.summarysize(stream) < 50_000
        finally
            close(stream)
        end
        sampled = tree_sampler(path;n=3,replace=false,rng=StableRNG(42))
        @test sampled.population_size == 5000 && length(sampled) == 3
        @test Base.summarysize(sampled) < 100_000
    end
end

@testset "Posterior trees: beta means, selected pairs and regionalization" begin
    stream_test_files() do dir,paths
        data = CommunityData(sparse([1 1 0; 1 0 1; 0 1 1; 0 0 0]);
            species_labels=["B","C","A"],sample_labels=["bc","ab","ac","empty"])
        trees = collect(read_trees(paths))
        for family in (:sorensen,:jaccard), component in (:total,:turnover,:nestedness)
            for tree in trees
                reverse_pair = getproperty(phylobeta(data,tree;index_family=family,pairs=[(1,2),(2,1)]),component)
                @test isequal(reverse_pair[1],reverse_pair[2])
            end
            expected = mean_dist((phylobeta(data,tree;index_family=family) for tree in trees);component)
            actual = read_trees(paths) do source
                mean_dist((phylobeta(data,tree;index_family=family) for tree in source);component)
            end
            @test actual.sample_labels == data.sample_labels && actual.replicates == 5
            @test actual.counts == expected.counts
            @test isapprox(actual.distances,expected.distances;nans=true)
            requested = [(1,2),(3,1),(2,1),(4,1),(4,4)]
            selected = read_trees(paths) do source
                mean_dist((phylobeta(data,tree;index_family=family,pairs=requested) for tree in source);
                    component,pairs=requested)
            end
            @test selected.counts == [actual.counts[i,j] for (i,j) in requested]
            @test isequal(selected.distances,[actual.distances[i,j] for (i,j) in requested])
        end
        defined = CommunityData(data.matrix[1:3,:];species_labels=data.species_labels,
            sample_labels=data.sample_labels[1:3])
        sample = tree_sampler(paths;n=8,rng=StableRNG(42))
        actual = mean_dist(phylobeta(defined,tree) for tree in sample)
        expected = mean_dist(phylobeta(defined,trees[i]) for i in sample.sample_indices)
        @test actual.counts == fill(8,3,3) && actual.distances ≈ expected.distances
        @test phyloregion(actual;k=2).membership == phyloregion(expected;k=2).membership
        @test prepare_phyloregion(actual).sample_labels == defined.sample_labels
    end
end
