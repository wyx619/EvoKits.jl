using Test, EvoKits, DataFrames, Random, SparseArrays, CSV

function extensions_tree(newick)
    mktempdir() do dir
        path = joinpath(dir, "tree.nwk")
        write_tree(path, newick)
        read_tree(path)
    end
end

@testset "Functional group selection and shared preparation" begin
    traits = DataFrame(species=string.('A':'F'), x=[0., 1, 2, 9, 10, 11], type=["a","a","a","b","b","b"])
    copy_traits = deepcopy(traits)
    rng = MersenneTwister(71)
    selection = optimal_functional_groups(traits; k=(4,2,3), nstart=12,
        numeric_scale=:none, gamma=1, rng)
    replay = MersenneTwister(71)
    fits = [functional_groups(traits; k, nstart=12, numeric_scale=:none, gamma=1, rng=replay) for k in 2:4]
    @test selection.k == [2,3,4]
    @test selection.objectives == [f.objective for f in fits]
    @test selection.iterations == [f.iterations for f in fits]
    @test selection.converged == [f.converged for f in fits]
    @test rand(rng) == rand(replay)
    @test selection.reason == :ok && selection.optimal_k == 3
    @test selection.grouping.membership == fits[2].membership
    @test isequal(traits, copy_traits)
    data = CommunityData(sparse([1 1 0 0 0 1; 0 0 1 1 1 0]); species_labels=traits.species)
    subset = optimal_functional_groups(data, traits[reverse(1:6), :]; k=2:4,
        numeric_scale=:none, gamma=1, nstart=12, rng=MersenneTwister(71))
    @test subset.objectives == selection.objectives
    @test functional_beta(data, subset.grouping).community.matrix ==
        collapse_traits(data, subset.grouping).matrix
    @test selection.grouping.gamma == 1
    @test selection.grouping.numeric_offsets == [0.]
    for (curve, reason) in (([0.,0.], :zero_objective),)
        index, why = EvoKits._functional_elbow(curve)
        @test index === nothing && why == reason
    end
    tied = DataFrame(species=["a","b","c"], x=[1,1,1])
    @test optimal_functional_groups(tied; k=1:3).reason == :zero_objective
    failed = optimal_functional_groups(traits; k=2:4, maxiter=1, rng=MersenneTwister(1))
    @test failed.reason == :nonconverged_candidates
    @test failed.optimal_k === nothing && failed.grouping === nothing
    partial = vcat(traits, DataFrame(species=["Z"], x=[missing], type=["z"]))
    @test optimal_functional_groups(partial; k=2:4, missing_traits=:skip).grouping.excluded_species == ["Z"]
    for bad in ((1,), (2,2), (0,2), (true,2), (2,100), 3, (2.0,3.0))
        @test_throws ArgumentError optimal_functional_groups(traits; k=bad)
    end
    for options in ((;low=0), (;low=NaN), (;max_fraction=1), (;nstart=0), (;maxiter=true))
        @test_throws ArgumentError optimal_functional_groups(traits; k=2:4, options...)
    end
end

function extensions_kl(p, q, method)
    p, q = BigFloat(p), BigFloat(q)
    poisson(a,b) = a == 0 ? b : b == 0 ? BigFloat(Inf) : a*log(a/b)+b-a
    value = poisson(p,q)
    method === :bernoulli && (value += poisson(1-p,1-q))
    Float64(value)
end

@testset "Indicator species KL formulas and ranking" begin
    theta = [0.9 0.1 0.2; 0.2 0.8 0.3; 0.1 0.3 0.7; 0.4 0.4 0.4]
    original = copy(theta)
    for method in (:poisson,:bernoulli), comparison in (:min,:max), shared in (false,true)
        result = indicators(theta; method, comparison, shared, top_indicators=5,
            species_labels=["a","b","c","d"], region_labels=["west","east","south"])
        @test size(result.indices) == (3,5)
        for region in 1:3
            reducer = comparison === :min ? minimum : maximum
            raw = [reducer([extensions_kl(theta[i,region], theta[i,j], method)
                for j in 1:3 if j != region]) for i in 1:4]
            eligible = [i for i in 1:4 if shared || argmax(theta[i,:]) == region]
            order = sort(eligible; by=i -> (-raw[i],i))
            @test result.counts[region] == length(order)
            @test result.indices[region,1:length(order)] == order
            @test result.scores[region,1:length(order)] ≈ raw[order]
            @test all(ismissing, result.indices[region,length(order)+1:end])
        end
    end
    @test theta == original
    for method in (:poisson,:bernoulli)
        limits = indicators([0.0 1; 1 0; 0 0; 1 1]; method, shared=true, top_indicators=4)
        @test all(x -> !ismissing(x) && !isnan(x), limits.scores)
        for region in 1:2, rank in 1:4
            i = limits.indices[region,rank]
            source, target = region == 1 ? ([0.,1,0,1], [1.,0,0,1]) : ([1.,0,0,1], [0.,1,0,1])
            @test limits.scores[region,rank] == extensions_kl(source[i], target[i], method)
        end
    end
    near = [0.5 nextfloat(0.5); 1e-200 2e-200]
    for method in (:poisson,:bernoulli)
        result = indicators(near; method, shared=true, top_indicators=2)
        for region in 1:2, rank in 1:2
            i = result.indices[region,rank]
            expected = extensions_kl(near[i,region],near[i,3-region],method)
            @test isapprox(result.scores[region,rank], expected; rtol=1e-12, atol=1e-290)
        end
    end
    @test indicators(ones(2,3); top_indicators=2).counts == [2,0,0]
    @test isequal(indicators(sparse(theta)).indices, indicators(theta).indices)
    for matrix in (zeros(0,2), zeros(2,1), [-1. 0; 0 1], [NaN 0; 0 1], [Inf 0; 0 1])
        @test_throws ArgumentError indicators(matrix)
    end
    @test_throws ArgumentError indicators([0. 2; 1 0]; method=:bernoulli)
    @test_throws ArgumentError indicators(theta; species_labels=fill("a",4))
    @test_throws ArgumentError indicators(theta; method=:normal)
    @test_throws ArgumentError indicators(theta; comparison=:mean)
    @test_throws ArgumentError indicators(theta; top_indicators=0)
end

function extensions_partition_sse(values, classes)
    sum(sum((BigFloat(x) - sum(BigFloat.(values[classes .== group])) /
        count(==(group),classes))^2 for x in values[classes .== group])
        for group in unique(classes))
end
function extensions_best_sse(values,k)
    sorted = sort(values)
    boundaries = findall(i -> sorted[i] != sorted[i+1], 1:length(sorted)-1)
    best = BigFloat(Inf)
    for mask in 0:(2^length(boundaries)-1)
        count_ones(mask) == k-1 || continue
        split = [boundaries[j] for j in eachindex(boundaries) if !iszero(mask & (1 << (j-1)))]
        classes = [searchsortedfirst(split,i) for i in eachindex(sorted)]
        best = min(best, extensions_partition_sse(sorted,classes))
    end
    best
end

@testset "Jenks exact partition, duplicates and stable variance" begin
    random = MersenneTwister(91)
    for trial in 1:18
        values = Float64.(rand(random,-5:15,8))
        for k in 1:min(4,length(unique(values)))
            result = classify_values(values; k)
            optimum = extensions_best_sse(values,k)
            actual = extensions_partition_sse(values,result.classes)
            @test isapprox(actual,optimum; rtol=1e-12,atol=1e-12)
            @test isapprox(result.objective,Float64(actual); rtol=1e-12,atol=1e-12)
            @test all(result.classes[i] == result.classes[j]
                for i in 1:8 for j in 1:8 if values[i] == values[j])
        end
    end
    shifted = 1e14 .+ [0.,0.03125,0.0625,10,10.03125,11]
    result = classify_values(shifted; k=2)
    @test result.classes == [1,1,1,2,2,2]
    @test result.objective ≈ Float64(extensions_partition_sse(shifted,result.classes))
    extremes = classify_values([-1e308,0.,1e308]; k=1)
    @test isinf(extremes.objective) && isfinite(extremes.normalized_objective)
    for method in (:jenks,:equal,:quantile)
        constant = classify_values(fill(3.,4); k=8, method)
        @test constant.actual_k == 1 && constant.reason == :reduced_classes
        @test constant.classes == ones(4)
        retained = classify_values([0.,NaN,1,missing,-2]; k=2, method, invalid=:skip)
        @test retained.classes[2] === missing && retained.classes[4] === missing
        @test retained.values[[1,3,5]] == [0,1,-2]
    end
    @test classify_values([0,1,2,3]; k=2,method=:equal).upper_bounds == [1.5,3.]
    @test classify_values([0,1,2,3]; k=2,method=:quantile).upper_bounds == [1.5,3.]
    for bad in ([],[NaN],[missing],[Inf])
        @test_throws ArgumentError classify_values(bad; invalid=:skip)
    end
    @test_throws ArgumentError classify_values(1:100; k=5,max_work=10)
    @test_throws ArgumentError classify_values([1,2]; k=0)
    @test_throws ArgumentError classify_values([1,2]; method=:foo)
end

@testset "Taxonomy-aided tree completion and provenance" begin
    tree = extensions_tree("((G_a:1,G_b:1)G:1,(H_a:1,I_a:1)Family:1)Root;")
    original = deepcopy(tree)
    table = DataFrame(species=["G_a","G_c","H_b","J_a","Z_a"],
        genus=["G","G","H","J","Z"], family=[missing,missing,missing,"Family",missing])
    built = phylobuilder(table,tree; extract=false)
    @test built.status == [:tree,:genus,:genus,:family,:missing]
    @test built.added_species == ["G_c","H_b","J_a"] && built.unmatched_species == ["Z_a"]
    @test built.tree.tip_labels == vcat(tree.tip_labels,built.added_species)
    @test all(==(2.),built.tree.dist_from_root[built.tree.tip_ids])
    @test pairwise_distance(built.tree)[1:4,1:4] == pairwise_distance(tree)
    @test built.tree.edge_length !== tree.edge_length && built.tree.children !== tree.children
    @test tree.tip_labels == original.tip_labels && tree.edge_length == original.edge_length
    extracted = phylobuilder(table,tree)
    @test Set(extracted.tree.tip_labels) == Set(["G_a","G_c","H_b","J_a"])
    @test extracted.status == built.status
    permuted = phylobuilder(table[reverse(1:5),:],tree; extract=false)
    for a in built.tree.tip_labels, b in built.tree.tip_labels
        @test distance(built.tree,a,b) ≈ distance(permuted.tree,a,b)
    end
    single = EvoKits._compact_tree_from_topology(Int32[0], [Int32[]], [0.],
        BitVector([true]), ["A_old"], ["A_old"])
    @test phylobuilder(["A_new"],single; extract=false).tree.ntips == 2
    nonmono = extensions_tree("((G_a:1,H_a:1):1,G_b:2);")
    @test_throws ArgumentError phylobuilder(["G_c"],nonmono)
    @test phylobuilder(["G_c"],nonmono; nonmonophyly=:mrca).status == [:genus]
    ambiguous = extensions_tree("((A_a:1,A_b:1)Family:1,(B_a:1,B_b:1)Family:1);")
    @test_throws ArgumentError phylobuilder(DataFrame(species=["Z_a"],family=["Family"]),ambiguous)
    uneven = extensions_tree("(G_a:1,G_b:2);")
    @test_throws ArgumentError phylobuilder(["G_c"],uneven)
    zero = phylobuilder(["G_c"],uneven; extract=false,branch_length=:zero)
    @test zero.tree.dist_from_root[zero.tree.tipname_to_id["G_c"]] == 0
    @test phylobuilder(["G_c"],uneven; extract=false,branch_length=3).tree.dist_from_root[end] == 3
    @test_throws ArgumentError phylobuilder(["G_a","G_a"],tree)
    @test_throws ArgumentError phylobuilder(["unknown"],tree)
    @test phylobuilder(["unknown"],tree; extract=false).status == [:missing]
    @test_throws ArgumentError phylobuilder(["G_c"],tree;branch_length=-1)
    @test_throws ArgumentError align_community(CommunityData([1;;];species_labels=["G_c"]),tree)
end

@testset "EDGE definitions and labeled categories" begin
    tree = extensions_tree("((A:1,B:1):1,(C:1,D:1):1);")
    result = edge(tree,["LC","NT","VU","CR"])
    @test result.ed == fill(1.5,4)
    @test result.ge == [0.001,0.01,0.1,0.999]
    @test result.edge ≈ log1p.(result.ed) .+ result.ge .* log(2)
    table = DataFrame(species=["D","A","C"],category=["CR","DD","EN"])
    subset = edge(tree,table)
    @test subset.status == [:unassessed,:absent_record,:ok,:ok]
    @test isequal(subset.edge[1:2],[missing,missing])
    custom = edge(tree,["LC","NT","VU","CR"];ge_mapping=Dict("LC"=>0,"NT"=>1,"VU"=>2,"CR"=>4))
    @test custom.edge ≈ log(2.5) .+ [0,1,2,4].*log(2)
    @test edge(tree,[missing,"DD","EW","EX"]).status == [:missing_category,:unassessed,:ok,:ok]
    @test_throws ArgumentError edge(tree,table;missing_category=:error)
    @test_throws ArgumentError edge(tree,["LC"];tip_labels=["Z"])
    @test_throws ArgumentError edge(tree,["LC","NT"];tip_labels=["A","A"])
    @test_throws ArgumentError edge(tree,["NE","NT","VU","CR"])
    @test_throws ArgumentError edge(tree,["LC","NT","VU","CR"];ge_mapping=Dict("LC"=>NaN))
end

@testset "Hotspots/coldspots and independent thresholds" begin
    values = [1.,2,3,4,5]
    @test hotspots(values;prob=20).threshold ≈ 4.2
    @test hotspots(values;prob=20).selected == [false,false,false,false,true]
    @test coldspots(values;prob=20).selected == [true,false,false,false,false]
    @test hotspots([1,2,2,2];prob=25).selected == [false,true,true,true]
    @test all(hotspots(fill(3.,5)).selected) && all(coldspots(fill(3.,5)).selected)
    result = hotspots([1.,NaN,missing,8]; reference=[0,4,8],prob=50,invalid=:skip)
    @test result.threshold == 4 && result.valid_count == 2 && result.reference_count == 3
    @test isequal(result.selected,[false,missing,missing,true])
    @test all(hotspots(values;prob=100).selected)
    @test hotspots(Float64[];reference=[1]).selected == Bool[]
    tree = extensions_tree("(A:1,B:1);")
    data = CommunityData([1 0; 1 1];species_labels=["A","B"],sample_labels=["n","s"])
    for metric in (pd(data,tree),weighted_endemism(data),phylo_endemism(data,tree))
        @test hotspots(metric).sample_labels == ["n","s"]
    end
    for prob in (0,-1,101,NaN,true)
        @test_throws ArgumentError hotspots(values;prob)
    end
    for value in ([],[NaN],[missing],[Inf])
        @test_throws ArgumentError coldspots(value;invalid=:skip)
    end
end

@testset "Compact local R reference" begin
    path = joinpath(@__DIR__,"fixtures","phyloregion_extensions","reference.tsv")
    if isfile(path)
        for row in CSV.File(path; delim='\t')
            values = parse.(Float64,split(row.input,','))
            if row.kind == "jenks"
                fit = classify_values(values;k=Int(row.option))
                upper = parse.(Float64,split(row.output,','))
                reference_classes = [searchsortedfirst(upper,x) for x in values]
                @test fit.objective ≈ Float64(extensions_partition_sse(values,reference_classes))
            elseif row.kind == "elbow"
                index, reason = EvoKits._functional_elbow(values)
                @test reason == :ok && index == row.option
            elseif row.kind == "hot"
                @test hotspots(values;prob=row.option).threshold ≈ parse(Float64,row.output)
            end
        end
    else
        error("Missing committed R reference fixture")
    end
end

