using EvoKits, Test, Random, StableRNGs, SparseArrays

function nexus_test_file(f)
    mktempdir() do dir
        path = joinpath(dir, "posterior.any")
        text = "\ufeff [preamble [nested; comment]]\n#nExUs\n" *
            "BEGIN TAXA;\nDIMENSIONS NTAX=3;\nTAXLABELS A B C;\nEND;\n" *
            "BeGiN TrEeS;\nTITLE 'posterior run';\nTRANSLATE\n  1 A,\n  2 B,\n  3 C;\n" *
            "TREE * STATE_0 = [&R] ((1:1[&rate=0.5],2:2)1:3,3:4);\n" *
            "tree 'STATE=10;run' = [&R] ((3:3,1:2):1,2:5);\n" *
            "UTREE STATE_20 = [&U] ((1:2,3:4):1,2:3);\nENDBLOCK;\n" *
            "BEGIN MRBAYES;\n[footer] set autoclose=yes;\nEND;\n"
        write(path, text)
        plain = joinpath(dir,"posterior.nwk")
        write(plain,"((A:1,B:2)1:3,C:4);\n((C:3,A:2):1,B:5);\n((A:2,C:4):1,B:3);")
        f(dir, path, plain)
    end
end

function nexus_error(path; kwargs...)
    stream = read_trees(path; kwargs...)
    error = try
        collect(stream)
        nothing
    catch caught
        caught
    end
    @test error isa ArgumentError
    @test !isopen(stream) && stream.io === nothing
    return sprint(showerror, error)
end

@testset "NEXUS: auto-detection, blocks, translation and filtering" begin
    nexus_test_file() do dir,path,plain
        stream = read_trees(path)
        @test stream.io === nothing
        trees, expected = collect(stream), collect(read_trees(plain))
        @test length(trees) == 3
        @test !isopen(stream) && stream.io === nothing
        for (actual,reference) in zip(trees,expected)
            @test actual.tip_labels == reference.tip_labels
            @test actual.node_labels == reference.node_labels
            @test actual.parent_of_edge == reference.parent_of_edge
            @test actual.child_of_edge == reference.child_of_edge
            @test actual.edge_length == reference.edge_length
            @test actual.tipname_to_id == reference.tipname_to_id
        end
        @test "1" in trees[1].node_labels # internal support/node labels are not translated
        @test distance(trees[1],"A","C") == 8
        @test only(collect(read_trees(path;burnin=1,thin=2))).tip_labels == ["C","A","B"]
        @test isempty(collect(read_trees(path;burnin=3)))
        combined = collect(read_trees([plain,path];burnin=1,thin=2))
        @test length(combined) == 2
        @test combined[1].tip_labels == combined[2].tip_labels == ["C","A","B"]
        @test distance(combined[1],"A","B") == distance(combined[2],"A","B")
        @test getproperty.(collect(read_trees([path,plain])),:tip_labels) ==
            getproperty.(vcat(trees,expected),:tip_labels)
        write(path,"#NEXUS\nBEGIN TREES; END;\nBEGIN TAXA; TAXLABELS A B; END;")
        @test isempty(collect(read_trees(path)))
        @test_throws ArgumentError tree_sampler(path;n=1)
        write(path,"#NEXUS [empty file]")
        @test isempty(collect(read_trees(path)))
        write(path,"#NEXUS[comment]BEGIN TREES;TREE t=(A:1,B:2);END;")
        @test only(collect(read_trees(path))).tip_labels == ["A","B"]
    end
end

@testset "NEXUS: quoted names, annotations and scoped tables" begin
    mktempdir() do dir
        path = joinpath(dir,"labels.nex")
        write(path,"#NEXUS\nBEGIN TREES; TRANSLATE 'code 1' 'A;[one]', " *
            "two 'B''quote', 三 '物种 C'; " *
            "TREE * 'tree=''one'';=' = [&R] (('code 1':1[&rate={1,2},note=\"a;b\"]," *
            "two:2)'two':3,三:4)[nested [x;y]]; END;")
        tree = only(collect(read_trees(path)))
        @test tree.tip_labels == ["A;[one]","B'quote","物种 C"]
        @test "two" in tree.node_labels
        @test distance(tree,"A;[one]","物种 C") == 8
        @test tree.node_labels[tree.tipname_to_id["B'quote"]] == "B'quote"
        # Tables reset on END, across TREES blocks, and across files.
        write(path,"#NEXUS BEGIN TREES; TRANSLATE 1 A,2 B;TREE t=(1:1,2:2);END;" *
            "BEGIN TREES;TRANSLATE 1 C,2 D;TREE u=(2:3,1:4);END;" *
            "BEGIN TREES;TREE v=(1:1,2:1);END;")
        trees = collect(read_trees(path))
        @test getproperty.(trees,:tip_labels) == [["A","B"],["D","C"],["1","2"]]
        kept = collect(read_trees(path;burnin=1,thin=2))
        @test length(kept) == 1 && kept[1].tip_labels == ["D","C"]
        other = joinpath(dir,"other.nex")
        write(other,"#NEXUS BEGIN TREES;TRANSLATE 1 E,2 F;TREE z=(1:5,2:6);END;")
        combined = collect(read_trees([path,other]))
        @test combined[4].tip_labels == ["E","F"]
        @test distance(combined[4],"E","F") == 11
        # Comments separate command/table tokens and can contain semicolons.
        write(path,"#NEXUS BEGIN[;]TREES; TRANSLATE 1[;]A,2[;]B; TREE[;]t[;]=[&R](1:1,2:2);END;")
        @test only(collect(read_trees(path))).tip_labels == ["A","B"]
        # Literal labels need no table, and underscore/spelling is preserved.
        write(path,"#NEXUS BEGIN TREES;TREE t=('G a':1,G_b:2);END;")
        @test only(collect(read_trees(path))).tip_labels == ["G a","G_b"]
        # Empty branch lengths retain the existing missing-length semantics.
        write(path,"#NEXUS BEGIN TREES;TRANSLATE 1 A,2 B;TREE t=(1,2);END;")
        @test all(isnan,only(collect(read_trees(path))).edge_length)
        # Commands/bodies can span lines; UTREE retains the written root.
        write(path,"#NEXUS\nBEGIN TREES;\nUTREE\n'unrooted example'\n= [&U]\n(A:1,B:2,C:3);\nEND;")
        sampled = tree_sampler(path;n=1)
        @test sampled.source_lines == [3] && sampled.source_indices == [1]
        tree = only(sampled)
        @test length(tree.children[tree.root]) == 3
        @test distance(tree,"A","C") == 4
        # Skipped non-tree blocks are not interpreted as posterior trees.
        write(path,"#NEXUS BEGIN DATA; DIMENSIONS NTAX=2 NCHAR=4; MATRIX 'A taxon' ACGA B ACCA; END;" *
            "BEGIN TREES;TREE t=('A taxon':1,B:2);END;")
        @test only(collect(read_trees(path))).tip_labels == ["A taxon","B"]
    end
end

@testset "NEXUS: errors, original tree positions and cleanup" begin
    mktempdir() do dir
        path = joinpath(dir,"bad.nex")
        for command in ("TRANSLATE;", "TRANSLATE 1 A,1 B;", "TRANSLATE 1 A,2 A;",
            "TRANSLATE 1 A,;", "TRANSLATE 1;", "TRANSLATE 1 A 2 B;", "TRANSLATE '' A;",
            "TRANSLATE 1 '';", "TRANSLATE 1 ' ';", "TRANSLATE '1'A;",
            "TRANSLATE 1 A;TRANSLATE 2 B;", "TREE t=(A:1,B:1);TRANSLATE 1 A;",
            "TREE t (A:1,B:1);", "TREE =(A:1,B:1);", "TREE * =(A:1,B:1);", "TREE t=;",
            "TREE ''=(A:1,B:1);", "TREE t extra=(A:1,B:1);", "LINK TAXA=other;",
            "UNSUPPORTED command;", "BEGIN OTHER;", "TITLE;", "'TREE' t=(A:1,B:1);" )
            write(path,"#NEXUS\nBEGIN TREES;\n" * command * "\nEND;")
            message = nexus_error(path;burnin=20)
            @test occursin(abspath(path),message)
            @test occursin("NEXUS file",message)
        end
        for text in ("#NEXUS\nBEGIN TREES;", "#NEXUS BEGIN TAXA; TAXLABELS A B;",
            "#NEXUS END;", "#NEXUS TREE t=(A:1,B:1);", "#NEXUS BEGIN TREES extra;",
            "#NEXUS BEGIN TREES;TREE t=(A:1,B:1);END extra;",
            "#NEXUS BEGIN TREES;TREE t=(A:1,B:1);[broken",
            "#NEXUS BEGIN TREES;TREE 'broken=(A:1,B:1);",
            "#NEXUS BEGIN TREES;TREE t=(A:1,B:1);END")
            write(path,text)
            nexus_error(path)
        end
        write(path,"#NEXUS\nBEGIN TAXA;\nTAXLABELS A B;\nEND;\nBEGIN TREES;\n" *
            "TRANSLATE 1 A,2 B;\nTREE first=(1:1,2:2);\nTREE second=(1:1,3:1);\nEND;")
        stream = read_trees(path)
        @test first(stream).tip_labels == ["A","B"]
        handle = stream.io
        message = try
            iterate(stream)
            ""
        catch error
            sprint(showerror,error)
        end
        @test occursin("tree 2, line 8",message)
        @test occursin("tip missing from TRANSLATE table: 3",message)
        @test !isopen(stream) && !isopen(handle)
        captured = Ref{TreeStream}()
        read_trees(path) do source
            captured[] = source
            first(source)
        end
        @test !isopen(captured[]) && isempty(captured[].nexus.translation)
        @test_throws ErrorException read_trees(path) do source
            first(source)
            error("stopped")
        end
        # A bad retained body fails lazily; burn-in may skip its structure.
        write(path,"#NEXUS\nBEGIN TREES;\nTREE bad=(A:1,,B:1);\nTREE good=(A:1,B:2);\nEND;")
        @test occursin("tree 1, line 3",nexus_error(path))
        @test only(collect(read_trees(path;burnin=1))).tip_labels == ["A","B"]
        @test only(tree_sampler(path;n=1,burnin=1)).ntips == 2
        @test_throws ArgumentError tree_sampler(path;n=2,replace=false)
        write(path,"#NEXUS BEGIN TREES;TRANSLATE 1 A,2 B;" *
            "TREE bad=(1:1,3:1);TREE good=(1:1,2:2);END;")
        seed = findfirst(seed -> rand(StableRNG(seed),1:2) == 2,1:20)
        @test tree_sampler(path;n=1,rng=StableRNG(seed)).source_indices == [2]
        @test only(collect(read_trees(path;burnin=1))).tip_labels == ["A","B"]
        @test_throws ArgumentError tree_sampler(path;n=2,replace=false)
        write(path,"#NEXUS BEGIN TREES;TREE t=(A:1,B:2)")
        @test_throws ArgumentError tree_sampler(path;n=1,burnin=10)
    end
end

@testset "NEXUS: sampling provenance, bounded state and posterior means" begin
    nexus_test_file() do dir,path,plain
        expected = collect(read_trees(plain))
        sampled = tree_sampler(path;n=40,rng=StableRNG(42))
        repeated = tree_sampler(path;n=40,rng=StableRNG(42))
        @test sampled.sample_indices == repeated.sample_indices
        @test sampled.population_size == 3
        @test Set(sampled.sample_indices) == Set(1:3)
        for draw in eachindex(sampled.trees)
            index = sampled.sample_indices[draw]
            actual, reference = sampled.trees[draw], expected[index]
            @test actual.tip_labels == reference.tip_labels
            @test actual.edge_length == reference.edge_length
            @test sampled.source_paths[draw] == abspath(path)
            @test sampled.source_indices[draw] == index
            @test sampled.source_lines[draw] == [13,14,15][index]
        end
        duplicate = findfirst(i -> i > 1 && sampled.sample_indices[i] == sampled.sample_indices[1],eachindex(sampled.trees))
        @test sampled.trees[1].tip_labels !== sampled.trees[duplicate].tip_labels
        @test sampled.trees[1].edge_length !== sampled.trees[duplicate].edge_length
        distinct = tree_sampler(path;n=3,replace=false,rng=StableRNG(7))
        @test sort(distinct.source_indices) == [1,2,3]
        @test_throws ArgumentError tree_sampler(path;n=4,replace=false)
        subset = tree_sampler([plain,path];n=2,replace=false,burnin=1,thin=2,rng=StableRNG(7))
        @test subset.population_size == 2 && subset.source_indices == [2,2]
        @test Set(subset.source_paths) == Set(abspath.([plain,path]))
        data = CommunityData(sparse([1 1 0;1 0 1;0 1 1;0 0 0]);
            species_labels=["B","C","A"],sample_labels=["bc","ab","ac","empty"])
        pairs = [(1,2),(2,1),(3,1),(4,4)]
        for family in (:sorensen,:jaccard), component in (:total,:turnover,:nestedness)
            average = read_trees(path) do trees
                mean_dist((phylobeta(data,tree;index_family=family) for tree in trees);component)
            end
            reference = read_trees(plain) do trees
                mean_dist((phylobeta(data,tree;index_family=family) for tree in trees);component)
            end
            @test isapprox(average.distances,reference.distances;nans=true)
            @test average.counts == reference.counts
            selected = read_trees(path) do trees
                mean_dist((phylobeta(data,tree;index_family=family,pairs) for tree in trees);component,pairs)
            end
            @test isequal(selected.distances,[average.distances[i,j] for (i,j) in pairs])
            @test selected.counts == [average.counts[i,j] for (i,j) in pairs]
        end
        defined = CommunityData(data.matrix[1:3,:];species_labels=data.species_labels,
            sample_labels=data.sample_labels[1:3])
        average = mean_dist(phylobeta(defined,tree) for tree in sampled)
        reference = mean_dist(phylobeta(defined,expected[i]) for i in sampled.sample_indices)
        @test average.distances ≈ reference.distances
        @test phyloregion(average;k=2).membership == phyloregion(reference;k=2).membership
        @test average.counts == fill(40,3,3)
        open(path,"w") do io
            println(io,"#NEXUS\nBEGIN TREES;TRANSLATE 1 A,2 B,3 C;")
            for i in 1:5000
                println(io,"TREE STATE_$i = [&R] ((1:1,2:2):3,3:4);")
            end
            println(io,"END;")
        end
        read_trees(path) do trees
            for _ in 1:2000
                iterate(trees)
            end
            @test Base.summarysize(trees) < 50_000
            @test length(trees.nexus.translation) == 3
        end
        result = tree_sampler(path;n=3,replace=false,rng=StableRNG(42))
        @test result.population_size == 5000
        @test Base.summarysize(result) < 100_000
        # A finished block releases a large translation table's capacity.
        open(path,"w") do io
            println(io,"#NEXUS\nBEGIN TREES;TRANSLATE")
            for i in 1:8000
                println(io,"$i unused_$i",i == 8000 ? ";" : ",")
            end
            println(io,"END;BEGIN TREES;TRANSLATE 1 A,2 B;TREE t=(1:1,2:2);END;")
        end
        read_trees(path) do trees
            @test first(trees).tip_labels == ["A","B"]
            @test Base.summarysize(trees) < 50_000
        end
    end
end
