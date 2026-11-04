using Test, EvoKits, LinearAlgebra, Statistics, StableRNGs, CSV, DataFrames

@testset "NMDS alignment: labeled rigid and similarity transforms" begin
    source = [0.0 0; 2 0; 0 3; 1 2; -2 1]
    labels = ["a", "b", "c", "d", "source_only"]
    q = [0.0 1; 1 0]
    reference = source[1:4, :] * q .+ [5.0 -3]
    reference = vcat(reference[[3,1,4,2], :], [20.0 10])
    target_labels = ["c", "a", "d", "b", "reference_only"]
    original = copy(source)
    fit = align_nmds(source, reference; sample_labels=labels, reference_labels=target_labels)
    @test fit isa NMDSAlignmentResult
    @test fit.sample_labels == labels && fit.sample_labels !== labels
    @test fit.matched_labels == labels[1:4]
    @test fit.source_indices == [1,2,3,4] && fit.reference_indices == [2,4,1,3]
    @test fit.points ≈ source * q .+ [5.0 -3] atol=2e-14
    @test fit.rotation ≈ q atol=1e-14
    @test fit.rotation' * fit.rotation ≈ Matrix{Float64}(I,2,2) atol=1e-14
    @test det(fit.rotation) ≈ -1
    @test fit.scale == 1 && !fit.scaled && fit.reflection
    @test fit.rmse < 1e-14 && maximum(fit.residuals) < 2e-14
    @test fit.source_center ≈ vec(mean(source[1:4,:]; dims=1))
    @test fit.points ≈ (source .- fit.source_center') * fit.rotation .* fit.scale .+ fit.reference_center'
    @test source == original
    for i in 1:5, j in 1:5
        @test norm(fit.points[i,:] - fit.points[j,:]) ≈ norm(source[i,:] - source[j,:]) atol=1e-14
    end
    proper = align_nmds(source, reference; sample_labels=labels,
        reference_labels=target_labels, reflection=false)
    @test det(proper.rotation) ≈ 1
    @test proper.rmse > 0.1 && proper.rmse >= fit.rmse
    scaled_reference = 2.7 .* reference
    scaled = align_nmds(source, scaled_reference; sample_labels=labels,
        reference_labels=target_labels, scale=true)
    @test scaled.scale ≈ 2.7
    @test scaled.points ≈ 2.7 .* (source * q .+ [5.0 -3]) atol=3e-14
    @test scaled.rmse < 3e-14 && scaled.scaled
    @test align_nmds(fit, reference; reference_labels=target_labels).points ≈ fit.points atol=1e-14
    for unit in (1e-200, 1e200)
        extreme = align_nmds(source .* unit, reference .* unit;
            sample_labels=labels, reference_labels=target_labels)
        @test extreme.rotation ≈ q atol=1e-14
        @test extreme.points ./ unit ≈ fit.points atol=2e-14
        @test isfinite(extreme.rmse)
    end
    line = reshape([-2.0, 0, 4], :, 1)
    one_axis = align_nmds(line, -line .+ 8)
    @test one_axis.points ≈ -line .+ 8
    @test one_axis.rotation ≈ fill(-1.0,1,1)
end

@testset "NMDS alignment: least-squares optimum and validation" begin
    rng = StableRNG(470)
    source = randn(rng, 9, 2)
    reference = randn(rng, 9, 2)
    fit = align_nmds(source, reference)
    x = source .- mean(source; dims=1)
    y = reference .- mean(reference; dims=1)
    # Independent planar formula; enumerate rotations and reflections.
    best = Inf
    for reflect in (1.0, -1.0)
        z = x * Diagonal([1.0, reflect])
        a = sum(z[:,1] .* y[:,1] + z[:,2] .* y[:,2])
        b = sum(z[:,2] .* y[:,1] - z[:,1] .* y[:,2])
        angle = atan(b, a)
        q = [cos(angle) -sin(angle); sin(angle) cos(angle)]
        best = min(best, sum(abs2, z*q-y))
    end
    @test fit.rmse^2 * 9 ≈ best atol=2e-13
    @test fit.rmse ≈ sqrt(mean(abs2, fit.residuals))
    @test fit.residuals ≈ [norm(row) for row in eachrow(fit.points-reference)]
    for bad in (zeros(0,2), zeros(3,0), fill(NaN,4,2), fill(Inf,4,2),
        fill(big"1e1000",4,2), fill(big"1e-1000",4,2), fill(missing,4,2))
        @test_throws ArgumentError align_nmds(bad, reference)
        @test_throws ArgumentError align_nmds(source, bad)
    end
    @test_throws ArgumentError align_nmds(source, reference[:,1:1])
    @test_throws ArgumentError align_nmds(source, reference; scale=1)
    @test_throws ArgumentError align_nmds(source, reference; reflection=1)
    @test_throws ArgumentError align_nmds(source, reference; sample_labels=fill("a",9))
    @test_throws ArgumentError align_nmds(source, reference; reference_labels=["a"])
    @test_throws ArgumentError align_nmds(source, reference; reference_labels=string.(11:19))
    @test_throws ArgumentError align_nmds(source[1:2,:], reference[1:2,:])
    @test_throws ArgumentError align_nmds(ones(4,2), source[1:4,:])
    @test_throws ArgumentError align_nmds(hcat(1:4, 1:4), source[1:4,:])
    @test_throws ArgumentError align_nmds(source[1:4,:], hcat(1:4, 1:4))
    @test_throws ArgumentError align_nmds(fit, reference; sample_labels=fit.sample_labels)
end

@testset "NMDS colors: palette, reusable bounds and constant axes" begin
    corners = [0.0 0; 0 1; 1 0; 1 1]
    labels = ["west", "east", "north", "south"]
    fit = nmds_colors(corners; sample_labels=labels)
    @test fit isa NMDSColorResult && fit.scale isa NMDSColorScale
    @test fit.hex == ["#DB9D85", "#86B875", "#CD99D8", "#4CB9CC"]
    @test fit.sample_labels == labels && fit.sample_labels !== labels
    @test fit.points == corners && fit.points !== corners
    @test fit.scale.lower == (0.0,0.0) && fit.scale.upper == (1.0,1.0)
    @test !any(fit.clipped)
    @test all(0 .<= fit.rgb .<= 1)
    for i in 1:4
        @test [parse(Int, fit.hex[i][a:(a+1)]; base=16) for a in (2,4,6)] ./ 255 == fit.rgb[i,:]
    end
    @test nmds_colors(corners .* 7 .- 20).hex == fit.hex
    shared = [0.1 0.2; 0.3 0.7; 0.8 0.6]
    reference = nmds_colors(vcat(corners, shared))
    subset = nmds_colors(shared[[3,1],:]; scale=reference.scale)
    @test subset.hex == reference.hex[[7,5]]
    @test !any(subset.clipped)
    @test nmds_colors(shared[[3,1],:]).hex != subset.hex
    outside = [-1.0 0.5; 0.5 2; 0.5 0.5]
    @test_throws ArgumentError nmds_colors(outside; scale=fit.scale)
    clipped = nmds_colors(outside; scale=fit.scale, outside=:clip)
    @test clipped.clipped == [true,true,false] && clipped.points == outside
    @test clipped.hex == nmds_colors(clamp.(outside,0,1); scale=fit.scale).hex
    constant = nmds_colors(fill(5.0,4,2))
    @test length(unique(constant.hex)) == 1
    @test constant.hex == nmds_colors(zeros(4,2)).hex
    @test nmds_colors(fill(0.0,1,2)).hex == constant.hex[1:1]
    @test nmds_colors(reshape([1.0,2,3],:,1)).hex == nmds_colors([1.0 0; 2 0; 3 0]).hex
    @test nmds_colors([1e300 0.0; 1e300 1e-300]).hex == nmds_colors([0.0 0; 0 1]).hex
    @test nmds_colors([-1e308 -1e308; 1e308 1e308]).hex == fit.hex[[1,4]]
    @test nmds_colors(corners .* 1e-300).hex == fit.hex
    @test corners == [0.0 0; 0 1; 1 0; 1 1]
end

@testset "NMDS postprocessing: result integration and invalid colors" begin
    positions = [0.0 0; 2 0; 0 1; 1 2; 3 3]
    distances = [norm(positions[i,:]-positions[j,:]) for i in 1:5, j in 1:5]
    labels = ["a", "b", "c", "d", "e"]
    base = nmds(distances; sample_labels=labels, nstart=1, init=positions, maxiter=0)
    before_points, before_stress = copy(base.points), base.stress
    fit = align_nmds(base, base)
    @test fit.points ≈ base.points atol=2e-15
    @test fit.sample_labels == labels
    colors = nmds_colors(base)
    @test colors.sample_labels == labels && colors.points == base.points
    @test base.points == before_points && base.stress == before_stress
    @test nmds_colors(fit; scale=colors.scale, outside=:clip).sample_labels == labels
    regions = phyloregion(distances; k=3, sample_labels=labels)
    regional = nmds(regions; nstart=1)
    regional_colors = nmds_colors(regional)
    @test regional_colors.sample_labels == string.(1:3)
    @test length(regional_colors.hex[regions.membership]) == length(labels)
    pair = nmds([0.0 1; 1 0]; dims=1, nstart=1)
    @test size(nmds_colors(pair).points) == (2,2)
    @test_throws ArgumentError nmds_colors(base; sample_labels=labels)
    @test_throws ArgumentError nmds_colors(base; outside=:ignore)
    @test_throws ArgumentError nmds_colors(base; scale=colors)
    @test_throws ArgumentError nmds_colors(zeros(4,3))
    @test_throws ArgumentError nmds_colors(zeros(0,2))
    @test_throws ArgumentError nmds_colors(fill(NaN,4,2))
    @test_throws ArgumentError nmds_colors(fill(Inf,4,2))
    @test_throws ArgumentError nmds_colors(fill(missing,4,2))
    @test_throws ArgumentError nmds_colors(positions; sample_labels=fill("a",5))
    @test_throws ArgumentError NMDSColorScale((1.0,0.0), (0.0,1.0))
    @test_throws ArgumentError NMDSColorScale((NaN,0.0), (1.0,1.0))
end

@testset "NMDS postprocessing: saved local R reference" begin
    table = CSV.read(joinpath(@__DIR__, "fixtures", "nmds_postprocess", "reference.tsv"),
        DataFrame; delim='\t', stringtype=String)
    unpack(text, n, d) = reshape(parse.(Float64, split(text, ',')), n, d)
    for row in eachrow(table)
        points = unpack(row.points, row.n, row.dims)
        labels = split(row.labels, '|')
        if row.kind == "alignment"
            reference = unpack(row.reference, row.reference_n, row.dims)
            fit = align_nmds(points, reference; sample_labels=labels,
                reference_labels=split(row.reference_labels, '|'), scale=row.scaled)
            @test fit.points ≈ unpack(row.expected, row.n, row.dims) atol=3e-13
            @test fit.rotation ≈ unpack(row.rotation, row.dims, row.dims) atol=3e-14
            @test fit.scale ≈ row.factor
            @test fit.rmse ≈ row.rmse atol=2e-14
        else
            fit = nmds_colors(points; sample_labels=labels)
            @test fit.hex == split(row.hex, '|')
            @test fit.sample_labels == labels
        end
    end
end
