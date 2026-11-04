using Test, EvoKits, SparseArrays, Random, CSV, DataFrames

_phylo_multi_same(x, y) = isequal(x, y) || isapprox(x, y; atol=5e-15, rtol=5e-13)

function phylo_multi_tree(newick)
    return mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, newick)
        read_tree(path)
    end
end

# Independent root walks, explicit branch set differences over every pair,
# and wider arithmetic. No production branch traversal, ranks or occupancy.
function phylo_multi_oracle(matrix, tree, family)
    return setprecision(BigFloat, 2304) do
        edge_of_child = Dict(tree.child_of_edge[e] => e for e in 1:tree.nedges)
        sets = [Set{Int}() for _ in axes(matrix, 1)]
        for site in axes(matrix, 1), tip in axes(matrix, 2)
            matrix[site, tip] > 0 || continue
            node = tree.tip_ids[tip]
            while node != tree.root
                push!(sets[site], edge_of_child[node])
                node = tree.parent_of_node[node]
            end
        end
        branch_length(edges) = sum(e -> BigFloat(tree.edge_length[e]), edges; init=BigFloat(0))
        values = branch_length.(sets)
        pooled = branch_length(union(sets...))
        a, m, x = sum(values) - pooled, BigFloat(0), BigFloat(0)
        for i in eachindex(sets), j in 1:(i-1)
            b, c = branch_length(setdiff(sets[i], sets[j])), branch_length(setdiff(sets[j], sets[i]))
            m += min(b, c)
            x += max(b, c)
        end
        denom = family === :sorensen ? 2a + m + x : a + m + x
        numerator = family === :sorensen ? m : 2m
        vdenom = a + numerator
        total = denom > 0 ? (m+x) / denom : BigFloat(NaN)
        turnover = vdenom > 0 ? numerator / vdenom : BigFloat(NaN)
        nestedness = total - turnover
        return (; total=Float64(total), turnover=Float64(turnover),
            nestedness=Float64(nestedness), pd=Float64.(values), pooled_pd=Float64(pooled))
    end
end

function check_phylo_multi(actual, expected; scalars_only=false)
    for field in (:total, :turnover, :nestedness, :pooled_pd)
        @test _phylo_multi_same(getproperty(actual, field), getproperty(expected, field))
    end
    if !scalars_only
        @test length(actual.pd) == length(expected.pd)
        @test all(_phylo_multi_same.(actual.pd, expected.pd))
    end
end

@testset "Phylogenetic multi-site beta: hand calculations and branch scope" begin
    tree = phylo_multi_tree("((A:1,B:2):3,(C:4,D:5):6);")
    matrix = [1 1 1 0; 1 0 0 1; 1 0 0 0]
    labels = ["ABC", "AD", "A"]
    for input in (matrix, Float64.(matrix), matrix .> 0, sparse(matrix),
        SparseMatrixCSC{Float64,Int32}(sparse(matrix)), view(matrix, :, :))
        sor = phylobeta_multi(input, tree; sample_labels=labels)
        jac = phylobeta_multi(input, tree; index_family=:jaccard)
        @test sor isa PhylobetaMultiResult
        @test sor.sample_labels == labels && sor.sample_labels !== labels
        @test sor.sample_indices == [1,2,3] && sor.richness == [3,2,1]
        @test sor.pd == [16,15,4] && sor.pooled_pd == 21
        @test sor.index_family === :sorensen && jac.index_family === :jaccard
        @test sor.total ≈ 17/31 && sor.turnover ≈ 5/19 && sor.nestedness ≈ 168/589
        @test jac.total ≈ 17/24 && jac.turnover ≈ 5/12 && jac.nestedness ≈ 7/24
        @test jac.sample_labels == ["sample_1", "sample_2", "sample_3"]
        for family in (:sorensen,:jaccard)
            selected = phylobeta_multi(input, tree; samples=[2,1], index_family=family)
            pair = phylobeta(input, tree; pairs=[(2,1)], index_family=family)
            for field in (:total,:turnover,:nestedness)
                @test _phylo_multi_same(getproperty(selected,field), only(getproperty(pair,field)))
            end
        end
    end
    # D and its whole ancestral clade are absent from this subset; keep the
    # original root-to-AB stem, but do not count any unoccupied CD branches.
    subset_matrix = [1 0 0 0; 0 1 0 0; 0 0 0 1]
    selected = phylobeta_multi(subset_matrix, tree; samples=[2,1], sample_labels=[11,22,33])
    @test selected.pd == [5,4] && selected.pooled_pd == 6
    @test selected.sample_indices == [2,1] && selected.sample_labels == ["22","11"]
    @test selected.total ≈ 1/3 && selected.turnover == 1/4 && selected.nestedness ≈ 1/12
    @test phylobeta_multi(subset_matrix, tree; samples=("22","11"),
        sample_labels=[11,22,33]).pd == selected.pd
    @test selected.pooled_pd < sum(tree.edge_length)
    @test subset_matrix == [1 0 0 0; 0 1 0 0; 0 0 0 1]
end

@testset "Phylogenetic multi-site beta: empty, zero, nested and disjoint sites" begin
    tree = phylo_multi_tree("((A:1,B:1):1,(C:1,D:1):1);")
    zero_tree = phylo_multi_tree("((A:0,B:0):0,(C:0,D:0):0);")
    star = phylo_multi_tree("((A:1,B:1):0,(C:1,D:1):0);")
    for family in (:sorensen,:jaccard)
        for (matrix, phy) in ((zeros(3,4),tree), ([1 0 0 0; 0 1 0 0],zero_tree))
            for input in (matrix,sparse(matrix))
                result = phylobeta_multi(input, phy; index_family=family)
                @test all(isnan,(result.total,result.turnover,result.nestedness))
                @test all(iszero,result.pd) && result.pooled_pd == 0
            end
        end
        one = phylobeta_multi([0 0 0 0; 1 0 0 0; 0 0 0 0],tree; index_family=family)
        @test one.total == 1 && isnan(one.turnover) && isnan(one.nestedness)
        mixed = phylobeta_multi([0 0 0 0; 1 0 0 0; 1 1 0 0],tree; index_family=family)
        @test mixed.turnover == 0 && mixed.total == mixed.nestedness
        identical = phylobeta_multi(ones(3,4),tree; index_family=family)
        @test identical.total == identical.turnover == identical.nestedness == 0
        nested = phylobeta_multi([1 1 1 0;1 1 0 0;1 0 0 0],tree; index_family=family)
        @test nested.turnover == 0 && nested.total == nested.nestedness
        matrix = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]
        disjoint = phylobeta_multi(matrix,star; index_family=family)
        @test disjoint.total == disjoint.turnover == 1 && disjoint.nestedness == 0
        taxonomic = beta_multi(matrix; index_family=family)
        @test all(f -> getproperty(disjoint,f) == getproperty(taxonomic,f),(:total,:turnover,:nestedness))
    end
end

@testset "Phylogenetic multi-site beta: independent oracle, order and context" begin
    rng = MersenneTwister(20261104)
    for trial in 1:30
        n, p = rand(rng,2:15), rand(rng,3:18)
        parts = ["t$i:$(rand(rng,0:12)/4)" for i in 1:p]
        while length(parts) > 2
            shuffle!(rng,parts)
            left, right = pop!(parts), pop!(parts)
            push!(parts,"($left,$right):$(rand(rng,0:12)/4)")
        end
        tree = phylo_multi_tree("($(join(parts,',')));")
        matrix = (rand(rng,n,p) .< rand(rng)) .* rand(rng,1:10,n,p)
        trial % 3 == 0 && (matrix[end,:] .= 0)
        trial % 4 == 0 && (matrix[1,:] .= matrix[2,:])
        indices = randperm(rng,n)[1:rand(rng,2:n)]
        perm = randperm(rng,p)
        labels = ["s$i" for i in 1:n]
        for family in (:sorensen,:jaccard), input in (matrix,sparse(matrix))
            context = prepare_community(input[:,perm],tree;tip_labels=tree.tip_labels[perm],sample_labels=labels)
            for subset in (collect(1:n),indices)
                expected = phylo_multi_oracle(matrix[subset,:],tree,family)
                result = phylobeta_multi(input[:,perm],tree;samples=subset,index_family=family,
                    tip_labels=tree.tip_labels[perm],sample_labels=labels)
                check_phylo_multi(result,expected)
                check_phylo_multi(phylobeta_multi(context;samples=subset,index_family=family),result)
                @test result.richness == vec(sum(matrix[subset,:] .> 0;dims=2))
                @test result.sample_indices == subset && result.sample_labels == labels[subset]
                @test all(v -> isnan(v) || 0 <= v <= 1,(result.total,result.turnover,result.nestedness))
                isfinite(result.nestedness) && @test result.total ≈ result.turnover + result.nestedness atol=5e-15
            end
            original = phylobeta_multi(context;samples=indices,index_family=family)
            reordered = phylobeta_multi(input[reverse(indices),:],tree;index_family=family)
            check_phylo_multi(reordered,original;scalars_only=true)
            @test reordered.pd ≈ reverse(original.pd)
        end
    end
end

@testset "Phylogenetic multi-site beta: long stems, close PD and extreme scales" begin
    matrix = [1 0 0; 0 1 0; 1 1 0; 0 0 0]
    for newick in ("((A:0.000001,B:0.000002):10000000000000000,C:1);",
        "((A:1e-308,B:2e-308):1e308,C:0);",
        "((A:5e-324,B:1e-323):1e308,C:0);",
        "((A:5e-324,B:1e-323):5e-324,C:0);",
        "((A:1,B:2):1e307,C:3);",
        "((A:0.1,B:0.3):0.2,C:0.4);")
        tree = phylo_multi_tree(newick)
        for rows in ([1,2],[1,2,3],[1,2,3,4]), family in (:sorensen,:jaccard)
            for input in (matrix,sparse(matrix))
                result = phylobeta_multi(input,tree;samples=rows,index_family=family)
                expected = phylo_multi_oracle(matrix[rows,:],tree,family)
                check_phylo_multi(result,expected)
                # Relative checks retain representable tiny components.
                for field in (:total,:turnover,:nestedness)
                    actual, target = getproperty(result,field), getproperty(expected,field)
                    isfinite(target) && @test actual ≈ target rtol=5e-13 atol=0
                end
            end
        end
        repeated = repeat(matrix[1:3,:],64,1)
        check_phylo_multi(phylobeta_multi(repeated,tree),phylo_multi_oracle(repeated,tree,:sorensen))
    end
end

@testset "Phylogenetic multi-site beta: ownership, alignment and invalid inputs" begin
    matrix = [1 1 0 0; 1 0 0 1; 0 1 0 0]
    for input in (matrix,sparse(matrix))
        tree = phylo_multi_tree("((A:1,B:2):3,(C:4,D:5):6);")
        labels, indices = ["north","south","middle"], [3,1]
        data = CommunityData(input[:,[4,1,2]];species_labels=["D","A","B"],sample_labels=labels)
        context = prepare_community(data,tree)
        result = phylobeta_multi(data,tree;samples=labels[indices])
        check_phylo_multi(result,phylobeta_multi(matrix,tree;samples=indices))
        @test result.sample_indices == indices && result.sample_labels == labels[indices]
        saved = deepcopy(context)
        from_context = phylobeta_multi(context;samples=indices)
        for field in (:pd,:richness,:sample_indices,:sample_labels)
            @test getproperty(from_context,field) !== getproperty(context,field === :pd ? :pd_values : field === :sample_indices ? :richness : field)
            fill!(getproperty(from_context,field),field === :sample_labels ? "changed" : 0)
        end
        @test all(field -> isequal(getproperty(context,field),getproperty(saved,field)),propertynames(context))
        fill!(tree.edge_length,Inf)
        fill!(tree.parent_of_node,0)
        fill!(data.matrix,0)
        check_phylo_multi(phylobeta_multi(context;samples=indices),result)
        @test labels == ["north","south","middle"] && indices == [3,1]
    end
    tree = phylo_multi_tree("((A:1,B:2):3,(C:4,D:5):6);")
    data = CommunityData(matrix;species_labels=tree.tip_labels)
    @test_throws ArgumentError phylobeta_multi(data,tree;sample_labels=nothing)
    @test_throws ArgumentError phylobeta_multi(data,tree;tip_labels=nothing)
    for selector in ([],[1],[1,1],[0,2],[1,4],[true,false],(1,2.0),(1,"2"),
        [missing,2],("missing","sample_1"),("sample_1","sample_1"),1,"sample_1")
        @test_throws ArgumentError phylobeta_multi(matrix,tree;samples=selector)
        @test_throws ArgumentError phylobeta_multi(prepare_community(matrix,tree);samples=selector)
    end
    for labels in (["a"],["a","a","c"],["a"," ","c"],[1,"1",3],[true,false,true])
        @test_throws ArgumentError phylobeta_multi(matrix,tree;sample_labels=labels)
    end
    for input in (zeros(0,4),zeros(1,4),spzeros(0,4),spzeros(1,4))
        @test_throws ArgumentError phylobeta_multi(input,tree)
    end
    @test_throws ArgumentError phylobeta_multi(matrix,tree;index_family=:wrong)
    @test_throws ArgumentError phylobeta_multi(prepare_community(matrix,tree);index_family=:wrong)
    @test_throws ArgumentError phylobeta_multi(ones(3,3),tree)
    @test_throws ArgumentError phylobeta_multi(matrix,tree;tip_labels=["A","A","C","D"])
    @test_throws ArgumentError phylobeta_multi(matrix,tree;tip_labels=["A","B","C","unknown"])
    @test_throws ArgumentError phylobeta_multi(fill(missing,2,4),tree)
    @test_throws ArgumentError phylobeta_multi(fill(1+im,2,4),tree)
    for bad in (-1.,Inf,-Inf,NaN), sparse_input in (false,true)
        input = Float64.(matrix)
        input[2,1] = bad
        @test_throws ArgumentError phylobeta_multi(sparse_input ? sparse(input) : input,tree;samples=[3,1])
        invalid_tree = deepcopy(tree)
        invalid_tree.edge_length[end] = bad
        @test_throws ArgumentError phylobeta_multi(matrix,invalid_tree;samples=[3,1])
    end
    overflow_tree = deepcopy(tree)
    fill!(overflow_tree.edge_length,1e308)
    @test_throws ArgumentError phylobeta_multi(matrix,overflow_tree)
    tiny = sparse([1,1,2,2,3],[1,2,1,3,2],[0.,nextfloat(0.),1e308,2.,1e-20],3,4)
    saved = copy(tiny)
    check_phylo_multi(phylobeta_multi(tiny,tree),phylobeta_multi(Matrix(tiny) .> 0,tree))
    @test tiny == saved && nonzeros(tiny) == nonzeros(saved) && nnz(tiny) == nnz(saved)
end

@testset "Phylogenetic multi-site beta: sparse and prepared storage" begin
    n,p = 6000,4000
    tree = phylo_multi_tree("($(join(("t$i:1" for i in 1:p),',')));")
    matrix = sparse(repeat(collect(1:n),3),repeat(1:3;inner=n),ones(Int,3n),n,p)
    phylobeta_multi(matrix,tree)
    bytes = @allocated result = phylobeta_multi(matrix,tree)
    @test result.total == result.turnover == result.nestedness == 0
    @test result.pooled_pd == 3 && all(==(3),result.pd)
    @test bytes < 20_000_000 # dense input: 192 MB; pair matrices: 864 MB
    @test Base.summarysize(result) < 1_000_000
    context = prepare_community(matrix,tree)
    phylobeta_multi(context)
    cached_bytes = @allocated cached = phylobeta_multi(context)
    check_phylo_multi(cached,result)
    @test cached_bytes < 5_000_000
end

@testset "Phylogenetic multi-site beta: saved local R reference" begin
    path = joinpath(@__DIR__,"fixtures","phylobeta_multi","reference.tsv")
    for row in eachrow(CSV.read(path,DataFrame;delim='\t',stringtype=String))
        matrix = reshape(parse.(Int,split(row.matrix,',')),row.n,row.p)
        columns = String.(split(row.columns,'|'))
        indices = parse.(Int,split(row.selected,'|'))
        tree = phylo_multi_tree(row.newick)
        result = phylobeta_multi(sparse(matrix),tree;tip_labels=columns,
            samples=indices,index_family=Symbol(row.family))
        check_phylo_multi(result,row;scalars_only=true)
    end
end
