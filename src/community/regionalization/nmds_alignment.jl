"""
    NMDSAlignmentResult

Coordinates returned by [`align_nmds`](@ref), in the source's original row
order. `matched_labels`, `source_indices` and `reference_indices` describe the
shared anchors in source order. The transform is
`scale * (source - source_center') * rotation + reference_center'`.
`residuals` are anchor distances and `rmse` is their root mean square, in
reference coordinate units. Arrays are owned; the original ordination and its
stress and convergence diagnostics are unchanged.
"""
struct NMDSAlignmentResult
    sample_labels::Vector{String}
    points::Matrix{Float64}
    matched_labels::Vector{String}
    source_indices::Vector{Int}
    reference_indices::Vector{Int}
    rotation::Matrix{Float64}
    scale::Float64
    source_center::Vector{Float64}
    reference_center::Vector{Float64}
    residuals::Vector{Float64}
    rmse::Float64
    scaled::Bool
    reflection::Bool
end

function _ordination_coordinates(input, labels, name)
    if input isa Union{NMDSResult, NMDSAlignmentResult}
        labels === nothing || throw(ArgumentError("$name already contains sample labels"))
        labels = input.sample_labels
        input = input.points
    end
    input isa AbstractMatrix{<:Real} ||
        throw(ArgumentError("$name must be an NMDS result or a real coordinate matrix without missing values"))
    Base.require_one_based_indexing(input)
    n, dims = size(input)
    n > 0 && dims > 0 || throw(ArgumentError("$name requires at least one row and coordinate axis"))
    points = Matrix{Float64}(input)
    all(isfinite, input) && all(isfinite, points) &&
        all(pair -> pair[1] == 0 || pair[2] != 0, zip(input, points)) ||
        throw(ArgumentError("$name coordinates must be finite Float64 values without underflow to zero"))
    return points, _validated_labels(labels, n, "sample")
end

function _ordination_center(points)
    unit = maximum(abs, points)
    unit > 0 || throw(ArgumentError("shared coordinates must have nonzero spread"))
    normalized = points ./ unit
    center = vec(mean(normalized; dims=1))
    normalized .-= transpose(center)
    return normalized, center, unit
end

function _ordination_full_rank(matrix)
    values = svdvals(matrix)
    return last(values) > max(size(matrix)...) * eps(Float64) * first(values)
end

"""
    align_nmds(source, reference; sample_labels=nothing, reference_labels=nothing,
               scale=false, reflection=true)

Align coordinates to a reference with orthogonal Procrustes least squares.
Inputs are `NMDSResult`, `NMDSAlignmentResult`, or real coordinate matrices.
Matrix labels default to `sample_1`, ...; result labels cannot be overridden.
Match exact shared labels and transform every source row, preserving its order.
Both inputs must have the same number `d` of axes, at least `d+1` shared labels,
and full-rank centered anchor configurations and cross-product. Insufficient
or degenerate anchors are errors; no arbitrary alignment is invented.

By default translation, rotation and reflection preserve all source distances.
`reflection=false` requires a proper rotation. `scale=true` also estimates one
positive uniform scale; it never stretches axes separately. Computation uses
scaled Float64 coordinates and rejects unrepresentable transforms or outputs.
Returns `NMDSAlignmentResult`, with the transform, anchor residuals and RMSE.
This aligns an existing configuration; it does not refit NMDS or change stress.
Shared labels must denote the same entities. Numeric region IDs from separate
clusterings are not automatically matched ecological identities.
"""
function align_nmds(source, reference; sample_labels=nothing, reference_labels=nothing,
    scale=false, reflection=true)
    scale isa Bool && reflection isa Bool ||
        throw(ArgumentError("scale and reflection must be Bool values"))
    points, labels = _ordination_coordinates(source, sample_labels, "source")
    target, target_labels = _ordination_coordinates(reference, reference_labels, "reference")
    dims = size(points, 2)
    size(target, 2) == dims || throw(ArgumentError("source and reference must have the same number of axes"))
    index = Dict(label => i for (i, label) in enumerate(target_labels))
    source_indices = [i for i in eachindex(labels) if haskey(index, labels[i])]
    reference_indices = [index[labels[i]] for i in source_indices]
    length(source_indices) > dims ||
        throw(ArgumentError("alignment requires at least d+1 shared labels for d coordinate axes"))
    left, source_center, source_unit = _ordination_center(points[source_indices, :])
    right, reference_center, reference_unit = _ordination_center(target[reference_indices, :])
    cross = transpose(left) * right
    all(_ordination_full_rank, (left, right, cross)) ||
        throw(ArgumentError("shared coordinates are rank deficient; select identifiable axes or additional anchors"))
    decomposition = svd(cross)
    signs = ones(dims)
    if !reflection && det(decomposition.U * decomposition.Vt) < 0
        signs[end] = -1
    end
    rotation = (decomposition.U .* transpose(signs)) * decomposition.Vt
    normalized_scale = scale ? dot(decomposition.S, signs) / sum(abs2, left) : 1.0
    multiplier = scale ? normalized_scale * reference_unit : source_unit
    fitted_scale = scale ? multiplier / source_unit : 1.0
    isfinite(fitted_scale) && fitted_scale > 0 && isfinite(multiplier) && multiplier > 0 ||
        throw(ArgumentError("alignment requires a positive representable Float64 scale; rescale coordinates"))
    aligned = (points ./ source_unit .- transpose(source_center)) * rotation
    aligned .*= multiplier
    reference_center .*= reference_unit
    source_center .*= source_unit
    aligned .+= transpose(reference_center)
    all(isfinite, aligned) || throw(ArgumentError("aligned coordinates exceed Float64 range; rescale coordinates"))
    differences = aligned[source_indices, :] - target[reference_indices, :]
    residuals = [norm(row) for row in eachrow(differences)]
    all(isfinite, residuals) || throw(ArgumentError("alignment residuals exceed Float64 range; rescale coordinates"))
    largest = maximum(residuals)
    rmse = largest == 0 ? 0.0 : (norm(residuals ./ largest) / sqrt(length(residuals))) * largest
    return NMDSAlignmentResult(labels, aligned, labels[source_indices], source_indices,
        reference_indices, rotation, fitted_scale, source_center, reference_center,
        residuals, rmse, scale, reflection)
end
