"""
    NMDSColorScale

Reusable two-axis coordinate bounds, available as `nmds_colors(...).scale`.
Both axes use the larger reference span, preserving their aspect ratio.
Constant reference axes map to the palette midpoint. Reuse this scale after
aligning subsequent configurations to the same reference.
"""
struct NMDSColorScale
    lower::NTuple{2,Float64}
    upper::NTuple{2,Float64}
    function NMDSColorScale(lower::NTuple{2,Float64}, upper::NTuple{2,Float64})
        all(isfinite, lower) && all(isfinite, upper) && all(lower .<= upper) ||
            throw(ArgumentError("color bounds must be finite and ordered on both axes"))
        return new(lower, upper)
    end
end

"""
    NMDSColorResult

Labeled output of [`nmds_colors`](@ref). `points` contains two coordinate axes
(one-axis input is padded with zero); `rgb` is quantized sRGB in [0,1] and
`hex` contains the matching uppercase `#RRGGBB` strings. `scale` retains the
reference bounds. `clipped` records rows outside those bounds when explicitly
allowed. All arrays are owned. Color differences are a visual summary of an
ordination, not a distance metric or classification confidence.
"""
struct NMDSColorResult
    sample_labels::Vector{String}
    points::Matrix{Float64}
    rgb::Matrix{Float64}
    hex::Vector{String}
    scale::NMDSColorScale
    clipped::BitVector
end

# Fixed L=70, C=50 polar CIELUV palette used by phyloregion::hexcols.
# D65 and the XYZ/sRGB conversion match the local R grDevices implementation.
# Only the two requested boundary colors are evaluated, not a 200×200 table.
function _nmds_hcl_color(hue)
    u = 50cosd(hue) / (13 * 70) + 0.1978398
    v = 50sind(hue) / (13 * 70) + 0.4683363
    y = 100 * ((70 + 16) / 116)^3
    x = 9y * u / (4v)
    z = -x / 3 - 5y + 3y / v
    linear = ((3.240479x - 1.537150y - 0.498535z) / 100,
        (-0.969256x + 1.875992y + 0.041556z) / 100,
        (0.055648x - 0.204043y + 1.057311z) / 100)
    return map(linear) do channel
        srgb = channel > 0.00304 ? 1.055 * channel^(1 / 2.4) - 0.055 : 12.92channel
        floor(Int, clamp(srgb, 0, 1) * 255 + 0.5)
    end
end

function _nmds_palette(x, y)
    # Preserve R's 200-level indexing and intermediate 8-bit boundary colors.
    a, b = ceil(Int, clamp(x, 0, 1) * 199) / 199, ceil(Int, clamp(y, 0, 1) * 199) / 199
    left = _nmds_hcl_color((1 - b) * 30 + b * 120)
    right = _nmds_hcl_color((1 - b) * 300 + b * 210)
    return ntuple(c -> floor(Int, (1 - a) * left[c] + a * right[c] + 0.5), 3)
end

function _nmds_color_position(value, axis, scale)
    lo, hi = scale.lower[axis], scale.upper[axis]
    lo == hi && return 0.5
    span = maximum(scale.upper .- scale.lower)
    if isfinite(span)
        return (value - lo) / span
    end
    unit = max(maximum(abs, scale.lower), maximum(abs, scale.upper))
    span = maximum(scale.upper ./ unit .- scale.lower ./ unit)
    return (value / unit - lo / unit) / span
end

"""
    nmds_colors(ordination; sample_labels=nothing, scale=nothing, outside=:error)

Map `NMDSResult`, `NMDSAlignmentResult`, or an n×1/n×2 real coordinate matrix
to labeled colors. One axis is padded with zero. Default bounds come from the
input, and both axes share the largest span. Constant axes use the palette
midpoint, including a single point or a completely constant configuration.

The fixed HCL boundary palette, intermediate sRGB rounding and 200-level
indexing follow local phyloregion `hexcols_rescale.R`. Only requested colors
are computed. No plotting, GIS or additional Julia dependency is required.
The mapping is sensitive to coordinate orientation. Use [`align_nmds`](@ref)
to compare fits with the same entities, then reuse the first result's `.scale`.
Scale reuse preserves the reference units and bounds; it does not align axes.
Out-of-bounds points are errors by default. `outside=:clip` explicitly clamps
their color coordinates and records them in `clipped`, preserving actual points.
Returns `NMDSColorResult`; inputs and reference scales are unchanged.
"""
function nmds_colors(ordination; sample_labels=nothing, scale=nothing, outside=:error)
    outside in (:error, :clip) || throw(ArgumentError("outside must be :error or :clip"))
    points, labels = _ordination_coordinates(ordination, sample_labels, "ordination")
    size(points, 2) in (1, 2) || throw(ArgumentError("color mapping requires one or two coordinate axes"))
    size(points, 2) == 1 && (points = hcat(points, zeros(size(points, 1))))
    reference = if scale === nothing
        NMDSColorScale(ntuple(a -> minimum(view(points, :, a)), 2),
            ntuple(a -> maximum(view(points, :, a)), 2))
    elseif scale isa NMDSColorScale
        scale
    else
        throw(ArgumentError("scale must be an NMDSColorScale from a reference color result"))
    end
    n = size(points, 1)
    rgb, hex, clipped = Matrix{Float64}(undef, n, 3), Vector{String}(undef, n), falses(n)
    for i in 1:n
        clipped[i] = any(a -> !(reference.lower[a] <= points[i, a] <= reference.upper[a]), 1:2)
        clipped[i] && outside === :error &&
            throw(ArgumentError("coordinates for $(labels[i]) exceed the reference color bounds; use outside=:clip explicitly"))
        positions = ntuple(a -> _nmds_color_position(
            clamp(points[i, a], reference.lower[a], reference.upper[a]), a, reference), 2)
        color = _nmds_palette(positions...)
        rgb[i, :] .= color ./ 255
        hex[i] = "#" * join(string(channel; base=16, pad=2) for channel in color) |> uppercase
    end
    return NMDSColorResult(labels, points, rgb, hex, reference, clipped)
end
