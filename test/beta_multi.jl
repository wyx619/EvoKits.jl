using Test, EvoKits, SparseArrays, Random, CSV, DataFrames

_multi_same(x, y) = isequal(x, y) || isapprox(x, y; atol=3e-15, rtol=3e-15)

# Independent set-difference enumeration and exact rational arithmetic.
# Nestedness is an exact difference here, unlike the production product.
function beta_multi_oracle(matrix, family)
    sets = [Set(findall(>(0), collect(row))) for row in eachrow(matrix)]
    gamma = length(union(sets...))
    a = BigInt(sum(length, sets) - gamma)
    m, x = big(0), big(0)
    for i in eachindex(sets), j in 1:(i-1)
        b, c = length(setdiff(sets[i], sets[j])), length(setdiff(sets[j], sets[i]))
        m += min(b, c)
        x += max(b, c)
    end
    d = family === :sorensen ? 2a + m + x : a + m + x
    v = family === :sorensen ? a + m : a + 2m
    total = d > 0 ? (m + x) // d : nothing
    turnover = v > 0 ? (family === :sorensen ? m : 2m) // v : nothing
    nestedness = total !== nothing && turnover !== nothing ? total - turnover : nothing
    converted(value) = value === nothing ? NaN : Float64(value)
    return (; total=converted(total), turnover=converted(turnover),
        nestedness=converted(nestedness), a, m, x, gamma, richness=length.(sets))
end

@testset "Multi-site beta: hand calculations and metadata" begin
    matrix = [1 1 1 0; 1 0 0 1; 1 0 0 0]
    labels = ["ABC", "AD", "A"]
    for input in (matrix, Float64.(matrix), matrix .> 0, sparse(matrix),
        SparseMatrixCSC{Float64,Int32}(sparse(matrix)), view(matrix, :, :))
        sor = beta_multi(input; sample_labels=labels)
        jac = beta_multi(input; index_family=:jaccard)
        @test sor isa BetaMultiResult
        @test sor.sample_labels == labels && sor.sample_labels !== labels
        @test sor.sample_indices == [1,2,3]
        @test sor.richness == [3,2,1] && sor.pooled_richness == 4
        @test sor.index_family === :sorensen && jac.index_family === :jaccard
        @test sor.total ≈ 3/5 && sor.turnover ≈ 1/3 && sor.nestedness ≈ 4/15
        @test jac.total ≈ 3/4 && jac.turnover ≈ 1/2 && jac.nestedness ≈ 1/4
        @test jac.sample_labels == ["sample_1", "sample_2", "sample_3"]
        for family in (:sorensen, :jaccard)
            selected = beta_multi(input; index_family=family, samples=[2,1])
            pair = beta_diss(input; index_family=family, pairs=[(2,1)])
            @test selected.sample_indices == [2,1]
            @test selected.richness == [2,3]
            for field in (:total, :turnover, :nestedness)
                @test _multi_same(getproperty(selected,field), only(getproperty(pair,field)))
            end
        end
    end
    # Selection owns its metadata, preserves request order and recomputes gamma.
    selected_rows = [3,1]
    fit = beta_multi(matrix; samples=selected_rows, sample_labels=labels)
    @test fit.sample_labels == ["A","ABC"] && fit.sample_indices == selected_rows
    @test fit.pooled_richness == 3 && fit.richness == [1,3]
    by_name = beta_multi(matrix; samples=("A","ABC"), sample_labels=labels)
    @test by_name.sample_indices == [3,1] && by_name.total == fit.total
    fit.sample_labels[1] = "changed"
    fit.sample_indices[1] = 2
    fit.richness[1] = 100
    @test labels == ["ABC","AD","A"] && selected_rows == [3,1]
    @test beta_multi(matrix; samples=selected_rows).richness == [1,3]
    @test matrix == [1 1 1 0; 1 0 0 1; 1 0 0 0]
end

@testset "Multi-site beta: empty, nested, disjoint and identical collections" begin
    for family in (:sorensen, :jaccard)
        for input in (zeros(3,4), spzeros(3,4), zeros(3,0), spzeros(3,0))
            result = beta_multi(input; index_family=family)
            @test all(isnan, (result.total,result.turnover,result.nestedness))
            @test result.pooled_richness == 0 && result.richness == [0,0,0]
        end
        one = beta_multi([0 0; 1 1; 0 0]; index_family=family)
        @test one.total == 1 && isnan(one.turnover) && isnan(one.nestedness)
        mixed = beta_multi([0 0; 1 0; 1 1]; index_family=family)
        @test mixed.turnover == 0 && mixed.total == mixed.nestedness
        @test mixed.total ≈ (family === :sorensen ? 2/3 : 4/5)
        nested = beta_multi([1 1 1; 1 1 0; 1 0 0]; index_family=family)
        @test nested.turnover == 0 && nested.total == nested.nestedness
        @test nested.total ≈ (family === :sorensen ? 2/5 : 4/7)
        disjoint = beta_multi([1 0 0; 0 1 0; 0 0 1; 0 0 0]; index_family=family)
        @test disjoint.total == disjoint.turnover == 1 && disjoint.nestedness == 0
        identical = beta_multi(ones(3,2); index_family=family)
        @test identical.total == identical.turnover == identical.nestedness == 0
    end
end

@testset "Multi-site beta: independent rational oracle and permutations" begin
    rng = MersenneTwister(20261104)
    for trial in 1:40
        n, p = rand(rng,2:20), rand(rng,0:35)
        matrix = (rand(rng,n,p) .< rand(rng)) .* rand(rng,1:9,n,p)
        trial % 3 == 0 && (matrix[end,:] .= 0)
        trial % 4 == 0 && (matrix[1,:] .= matrix[2,:])
        indices = randperm(rng,n)[1:rand(rng,2:n)]
        colperm = randperm(rng,p)
        for family in (:sorensen,:jaccard), input in (matrix,sparse(matrix))
            for subset in (collect(1:n),indices)
                oracle = beta_multi_oracle(matrix[subset,:],family)
                actual = beta_multi(input; samples=subset,index_family=family)
                @test actual.richness == oracle.richness && actual.pooled_richness == oracle.gamma
                for field in (:total,:turnover,:nestedness)
                    @test _multi_same(getproperty(actual,field),getproperty(oracle,field))
                end
                @test all(v -> isnan(v) || 0 <= v <= 1,
                    (actual.total,actual.turnover,actual.nestedness))
                if isfinite(actual.nestedness)
                    @test actual.total ≈ actual.turnover + actual.nestedness atol=3e-15
                end
            end
            reordered = beta_multi(input[reverse(indices),colperm]; index_family=family)
            original = beta_multi(input; samples=indices,index_family=family)
            @test all(field -> _multi_same(getproperty(reordered,field),getproperty(original,field)),
                (:total,:turnover,:nestedness))
        end
        oracle = beta_multi_oracle(matrix,:sorensen)
        stats = EvoKits._beta_multi_statistics(vec(sum(matrix .> 0;dims=2)),
            vec(sum(matrix .> 0;dims=1)))
        @test (stats.shared_excess,stats.minimum_unique,stats.maximum_unique,stats.pooled) ==
            (oracle.a,oracle.m,oracle.x,oracle.gamma)
    end
end

@testset "Multi-site beta: labels, CommunityData and input validation" begin
    matrix = [1 1 1 0; 1 0 0 1; 1 0 0 0]
    for input in (matrix,sparse(matrix))
        data = CommunityData(input; sample_labels=[11,22,33])
        saved = copy(data.matrix)
        fit = beta_multi(data; samples=["33","11"])
        @test fit.sample_labels == ["33","11"] && fit.sample_indices == [3,1]
        @test fit.richness == [1,3] && fit.pooled_richness == 3
        @test fit.total == beta_multi(matrix; samples=[3,1]).total
        @test data.matrix == saved && data.sample_labels == ["11","22","33"]
        @test_throws ArgumentError beta_multi(data; sample_labels=nothing)
        @test_throws ArgumentError beta_multi(data; tip_labels=["a","b","c","d"])
        for selector in ([],[1],[1,1],[0,2],[1,4],[true,false],(1,2.0),
            (1,"2"),[missing,2],("missing","sample_1"),("sample_1","sample_1"),1,"sample_1")
            @test_throws ArgumentError beta_multi(input; samples=selector)
        end
        @test_throws ArgumentError beta_multi(input; index_family=:wrong)
        for labels in (["a"],["a","a","c"],["a"," ","c"],[1,"1",3],[true,false,true])
            @test_throws ArgumentError beta_multi(input; sample_labels=labels)
        end
    end
    for input in (zeros(0,3),zeros(1,3),spzeros(0,3),spzeros(1,3))
        @test_throws ArgumentError beta_multi(input)
    end
    for invalid in (-1.0,NaN,Inf,-Inf), sparse_input in (false,true)
        input = Float64.(matrix)
        input[2,1] = invalid
        input = sparse_input ? sparse(input) : input
        @test_throws ArgumentError beta_multi(input; samples=[3,1]) # excluded rows still validated
    end
    @test_throws ArgumentError beta_multi(fill(missing,2,2))
    @test_throws ArgumentError beta_multi(fill(1+im,2,2))
    tiny = sparse([1,1,2,2,3],[1,2,1,3,2],
        [0.0,nextfloat(0.0),1e308,2.0,1e-20],3,4)
    saved = copy(tiny)
    expected = beta_multi(Matrix(tiny) .> 0)
    actual = beta_multi(tiny)
    @test actual.richness == [1,2,1] && actual.pooled_richness == 3
    @test all(f -> _multi_same(getproperty(actual,f),getproperty(expected,f)),(:total,:turnover,:nestedness))
    @test tiny == saved && nonzeros(tiny) == nonzeros(saved) && nnz(tiny) == nnz(saved)
    for input in ([big"1e-1000" 0; big"1e1000" 1], [1//2 0; 1//3 1])
        @test beta_multi(input).total == beta_multi([1 0; 1 1]).total
    end
end

@testset "Multi-site beta: integer cancellation and sparse storage" begin
    stats = EvoKits._beta_multi_statistics(fill(10_000,1000),fill(1000,10_000))
    @test stats.shared_excess == 9_990_000 && stats.minimum_unique == stats.maximum_unique == 0
    for family in (:sorensen,:jaccard), (a,m,x) in
        ((Int128(2)^80,Int128(2)^80,Int128(2)^80+1),(Int128(1),Int128(0),Int128(2)^80))
        total,turnover,nestedness = EvoKits._beta_multi_components(a,m,x,family)
        d = family === :sorensen ? 2big(a)+m+x : big(a)+m+x
        v = family === :sorensen ? big(a)+m : big(a)+2m
        exact_total = (big(m)+x)//d
        exact_turnover = (family === :sorensen ? big(m) : 2big(m))//v
        @test total ≈ Float64(exact_total)
        @test turnover ≈ Float64(exact_turnover)
        @test nestedness ≈ Float64(exact_total-exact_turnover)
        @test nestedness > 0
    end
    n,p = 6000,40_000
    matrix = sparse(repeat(collect(1:n),3),repeat(1:3;inner=n),ones(Int,3n),n,p)
    beta_multi(matrix)
    bytes = @allocated result = beta_multi(matrix)
    @test result.total == result.turnover == result.nestedness == 0
    @test bytes < 8_000_000 # dense input: 1.92 GB; three pair matrices: 864 MB
    @test Base.summarysize(result) < 1_000_000
end

@testset "Multi-site beta: saved local R reference" begin
    path = joinpath(@__DIR__,"fixtures","beta_multi","reference.tsv")
    for row in eachrow(CSV.read(path,DataFrame;delim='\t',stringtype=String))
        matrix = reshape(parse.(Int,split(row.matrix,',')),row.n,row.p)
        result = beta_multi(sparse(matrix); index_family=Symbol(row.family))
        for field in (:total,:turnover,:nestedness)
            @test _multi_same(getproperty(result,field),getproperty(row,field))
        end
        @test result.pooled_richness == row.pooled_richness
    end
end
