"""Reusable species-to-functional-group mapping returned by [`functional_groups`](@ref)."""
struct FunctionalGroupingResult
    species_labels::Vector{String}
    membership::Vector{Int}
    group_labels::Vector{String}
    group_sizes::Vector{Int}
    numeric_traits::Vector{String}
    categorical_traits::Vector{String}
    numeric_offsets::Vector{Float64}
    numeric_scales::Vector{Float64}
    numeric_centers::Matrix{Float64}
    categorical_modes::Matrix{String}
    numeric_weights::Vector{Float64}
    categorical_weights::Vector{Float64}
    numeric_scale::Symbol
    gamma::Float64
    objective::Float64
    iterations::Int
    converged::Bool
    nstart::Int
    excluded_species::Vector{String}
end

"""Functional beta, the aggregated community, and its owned species grouping."""
struct FunctionalBetaResult{B<:BetaDissResult}
    community::CommunityData
    grouping::FunctionalGroupingResult
    beta::B
    statistic::Symbol
    excluded_species::Vector{String}
end

function _functional_columns(table, species_col, trait_cols, categorical_traits)
    (species_col isa Symbol || species_col isa AbstractString) ||
        throw(ArgumentError("species_col must name a table column"))
    names_table = String.(names(table))
    species = String(species_col)
    species in names_table || throw(ArgumentError("traits table is missing species column $species"))
    columns = if trait_cols === nothing
        filter(!=(species), names_table)
    elseif trait_cols isa Symbol || trait_cols isa AbstractString
        [String(trait_cols)]
    elseif trait_cols isa AbstractVector || trait_cols isa Tuple
        all(x -> x isa Symbol || x isa AbstractString, trait_cols) ||
            throw(ArgumentError("trait_cols must contain column names"))
        String.(collect(trait_cols))
    else
        throw(ArgumentError("trait_cols must contain column names"))
    end
    !isempty(columns) && allunique(columns) && species ∉ columns &&
        all(in(names_table), columns) ||
        throw(ArgumentError("trait_cols must select distinct existing trait columns"))
    forced = if categorical_traits isa Symbol || categorical_traits isa AbstractString
        [String(categorical_traits)]
    elseif categorical_traits isa AbstractVector || categorical_traits isa Tuple
        all(x -> x isa Symbol || x isa AbstractString, categorical_traits) ||
            throw(ArgumentError("categorical_traits must contain column names"))
        String.(collect(categorical_traits))
    else
        throw(ArgumentError("categorical_traits must contain column names"))
    end
    allunique(forced) && all(in(columns), forced) ||
        throw(ArgumentError("categorical_traits must be distinct selected traits"))
    return species, columns, Set(forced)
end

function _functional_category(value, column)
    if value isa AbstractString || value isa Symbol
        label = String(value)
        isempty(strip(label)) && throw(ArgumentError("categorical trait $column cannot be empty"))
        return label, :text
    elseif value isa Bool
        return string(value), :bool
    elseif value isa Real
        isfinite(value) || throw(ArgumentError("categorical trait $column must be finite"))
        return string(value), (:number, typeof(value))
    end
    throw(ArgumentError("categorical trait $column must contain strings, symbols, booleans or finite numbers"))
end

function _functional_prepare(table::AbstractDataFrame, selected;
    species_col, trait_cols, categorical_traits, missing_traits,
)
    _community_missing_mapping(missing_traits, "missing_traits")
    species_col, columns, forced = _functional_columns(table, species_col, trait_cols, categorical_traits)
    labels = _validated_labels(table[!, species_col], nrow(table), "species")
    positions = Dict(label => i for (i, label) in enumerate(labels))
    requested = selected === nothing ? labels : selected
    rows, species, excluded = Int[], String[], String[]
    for label in requested
        row = get(positions, label, 0)
        missing_row = row == 0 || any(ismissing(table[row, column]) for column in columns)
        if missing_row
            missing_traits === :error && throw(ArgumentError("missing traits for species: $label"))
            push!(excluded, label)
        else
            push!(rows, row)
            push!(species, label)
        end
    end
    isempty(rows) && throw(ArgumentError("no species have complete selected traits"))
    numeric_columns, categorical_columns = String[], String[]
    for column in columns
        values = table[rows, column]
        if column ∉ forced && all(x -> x isa Real && !(x isa Bool), values)
            push!(numeric_columns, column)
        else
            push!(categorical_columns, column)
        end
    end
    numeric = Matrix{Float64}(undef, length(rows), length(numeric_columns))
    categorical = Matrix{String}(undef, length(rows), length(categorical_columns))
    for (j, column) in enumerate(numeric_columns), (i, row) in enumerate(rows)
        value = table[row, column]
        isfinite(value) || throw(ArgumentError("numeric trait $column must be finite"))
        converted = Float64(value)
        isfinite(converted) && (iszero(value) || !iszero(converted)) ||
            throw(ArgumentError("numeric trait $column must be representable as finite Float64"))
        numeric[i, j] = converted
    end
    for (j, column) in enumerate(categorical_columns)
        kind = nothing
        for (i, row) in enumerate(rows)
            category, current = _functional_category(table[row, column], column)
            (kind === nothing || kind === current) ||
                throw(ArgumentError("categorical trait $column mixes incompatible value types"))
            kind = current
            categorical[i, j] = category
        end
    end
    return species, excluded, numeric, categorical, numeric_columns, categorical_columns
end

function _functional_weights(weights, numeric_columns, categorical_columns)
    columns = vcat(numeric_columns, categorical_columns)
    weights === nothing && return ones(length(numeric_columns)), ones(length(categorical_columns))
    weights isa AbstractDict || throw(ArgumentError("trait_weights must be a dictionary keyed by trait name"))
    assigned = Dict{String, Float64}()
    for (key, value) in weights
        (key isa Symbol || key isa AbstractString) || throw(ArgumentError("trait weight names must be strings or symbols"))
        name = String(key)
        name in columns && !haskey(assigned, name) || throw(ArgumentError("unknown or repeated trait weight: $name"))
        value isa Real && !(value isa Bool) && isfinite(value) && value > 0 ||
            throw(ArgumentError("trait weights must be finite positive real numbers"))
        converted = Float64(value)
        isfinite(converted) && converted > 0 || throw(ArgumentError("trait weight is not representable as Float64"))
        assigned[name] = converted
    end
    return [get(assigned, c, 1.0) for c in numeric_columns],
        [get(assigned, c, 1.0) for c in categorical_columns]
end

function _functional_transform(raw, scale)
    scale in (:none, :sd, :range) || throw(ArgumentError("numeric_scale must be :none, :sd or :range"))
    n, q = size(raw)
    offsets, scales = zeros(q), ones(q)
    numeric = similar(raw)
    for j in 1:q
        offsets[j] = raw[1, j]
        delta = [raw[i, j] - offsets[j] for i in 1:n]
        all(isfinite, delta) || throw(ArgumentError("numeric trait differences overflowed; rescale the data"))
        if scale === :sd && n > 1
            avg = sum(x / n for x in delta)
            residual = delta .- avg
            all(isfinite, residual) || throw(ArgumentError("numeric trait variation overflowed; rescale the data"))
            largest = maximum(abs, residual)
            if largest > 0
                scales[j] = largest * sqrt(sum((x / largest)^2 for x in residual) / (n - 1))
            end
        elseif scale === :range
            scales[j] = maximum(delta) - minimum(delta)
            iszero(scales[j]) && (scales[j] = 1.0)
        end
        isfinite(scales[j]) && scales[j] > 0 ||
            throw(ArgumentError("numeric trait scale is not finite and positive"))
        for i in 1:n
            numeric[i, j] = delta[i] / scales[j]
            isfinite(numeric[i, j]) || throw(ArgumentError("scaled trait overflowed"))
        end
    end
    return numeric, offsets, scales
end

function _functional_gamma(numeric, categorical, nw, cw, gamma)
    if gamma !== :auto
        gamma isa Real && !(gamma isa Bool) && isfinite(gamma) && gamma > 0 ||
            throw(ArgumentError("gamma must be a finite positive real number or :auto"))
        result = Float64(gamma)
        isfinite(result) && result > 0 || throw(ArgumentError("gamma is not representable as Float64"))
        return result
    end
    n = size(numeric, 1)
    (n < 2 || isempty(nw) || isempty(cw)) && return 1.0
    numeric_pairs = 0.0
    for j in eachindex(nw)
        avg = sum(numeric[i, j] / n for i in 1:n)
        numeric_pairs += nw[j] * n * sum((numeric[i, j] - avg)^2 for i in 1:n)
    end
    categorical_pairs = 0.0
    for j in eachindex(cw)
        counts = Dict{String, Int}()
        for i in 1:n
            value = categorical[i, j]
            counts[value] = get(counts, value, 0) + 1
        end
        categorical_pairs += cw[j] * (Float64(n)^2 - sum(Float64(v)^2 for v in values(counts))) / 2
    end
    isfinite(numeric_pairs) && isfinite(categorical_pairs) ||
        throw(ArgumentError("automatic gamma calculation overflowed; rescale traits or set gamma"))
    (numeric_pairs > 0 && categorical_pairs > 0) || return 1.0
    result = numeric_pairs / categorical_pairs
    isfinite(result) && result > 0 || throw(ArgumentError("automatic gamma is not representable; set gamma"))
    return result
end

@inline function _functional_distance(numeric, categorical, row, means, modes, group, nw, penalties)
    distance = 0.0
    for j in eachindex(nw)
        difference = (numeric[row, j] - means[group, j]) * sqrt(nw[j])
        distance += difference * difference
    end
    for j in eachindex(penalties)
        categorical[row, j] == modes[group, j] || (distance += penalties[j])
    end
    isfinite(distance) || throw(ArgumentError("functional trait distance overflowed; rescale traits or weights"))
    return distance
end

function _functional_run(numeric, categorical, k, nw, penalties, rng, maxiter)
    n, q = size(numeric)
    r = size(categorical, 2)
    initial = randperm(rng, n)[1:k]
    means = copy(numeric[initial, :])
    modes = copy(categorical[initial, :])
    membership, previous, counts = zeros(Int, n), zeros(Int, n), zeros(Int, k)
    distances = zeros(n)
    iterations, converged = 0, false
    for iteration in 1:maxiter
        iterations = iteration
        copyto!(previous, membership)
        fill!(counts, 0)
        for row in 1:n
            best = previous[row] == 0 ? 1 : previous[row]
            best_distance = _functional_distance(numeric, categorical, row, means, modes, best, nw, penalties)
            for group in 1:k
                group == best && continue
                distance = _functional_distance(numeric, categorical, row, means, modes, group, nw, penalties)
                if distance < best_distance
                    best, best_distance = group, distance
                end
            end
            membership[row], distances[row] = best, best_distance
            counts[best] += 1
        end
        for empty_group in 1:k
            counts[empty_group] == 0 || continue
            donor, farthest = 0, -1.0
            for row in 1:n
                if counts[membership[row]] > 1 && distances[row] > farthest
                    donor, farthest = row, distances[row]
                end
            end
            donor > 0 || error("cannot repair an empty functional group")
            counts[membership[donor]] -= 1
            membership[donor] = empty_group
            counts[empty_group] += 1
            distances[donor] = 0.0
        end
        fill!(means, 0.0)
        for row in 1:n, j in 1:q
            group = membership[row]
            means[group, j] += numeric[row, j] / counts[group]
        end
        for group in 1:k, j in 1:r
            frequencies = Dict{String, Int}()
            for row in 1:n
                membership[row] == group || continue
                value = categorical[row, j]
                frequencies[value] = get(frequencies, value, 0) + 1
            end
            best_value, best_count = "", -1
            for (value, count) in frequencies
                if count > best_count || (count == best_count && isless(value, best_value))
                    best_value, best_count = value, count
                end
            end
            modes[group, j] = best_value
        end
        if membership == previous
            converged = true
            break
        end
    end
    objective = sum(_functional_distance(numeric, categorical, i, means, modes, membership[i], nw, penalties)
        for i in 1:n)
    isfinite(objective) || throw(ArgumentError("functional clustering objective overflowed"))
    return membership, means, modes, objective, iterations, converged
end

function _functional_workspace(table, selected; species_col, trait_cols, categorical_traits,
    missing_traits, numeric_scale, trait_weights, gamma,
)
    species, excluded, raw, categorical, nc, cc = _functional_prepare(table, selected;
        species_col, trait_cols, categorical_traits, missing_traits)
    numeric, offsets, scales = _functional_transform(raw, numeric_scale)
    nw, cw = _functional_weights(trait_weights, nc, cc)
    resolved_gamma = _functional_gamma(numeric, categorical, nw, cw, gamma)
    penalties = resolved_gamma .* cw
    all(isfinite, penalties) && all(>(0), penalties) || throw(ArgumentError("categorical penalties overflowed or underflowed"))
    return (; species, excluded, numeric, categorical, nc, cc, offsets, scales, nw, cw,
        numeric_scale, resolved_gamma, penalties)
end

function _functional_fit_prepared(work; k, nstart, maxiter, rng)
    (; species, excluded, numeric, categorical, nc, cc, offsets, scales, nw, cw,
        numeric_scale, resolved_gamma, penalties) = work
    k <= length(species) || throw(ArgumentError("k cannot exceed the number of species with complete traits"))
    best = nothing
    for _ in 1:Int(nstart)
        candidate = _functional_run(numeric, categorical, Int(k), nw, penalties, rng, Int(maxiter))
        if best === nothing || candidate[4] < best[4] ||
            (candidate[4] == best[4] && candidate[6] && !best[6])
            best = candidate
        end
    end
    membership, centers, modes, objective, iterations, converged = best
    # Canonicalize labels by the first species in each group.
    order = unique(membership)
    remap = zeros(Int, Int(k))
    for (new, old) in enumerate(order)
        remap[old] = new
    end
    membership = [remap[group] for group in membership]
    sizes = [count(==(group), membership) for group in 1:Int(k)]
    return FunctionalGroupingResult(copy(species), membership, ["functional_$i" for i in 1:Int(k)],
        sizes, copy(nc), copy(cc), copy(offsets), copy(scales), centers[order, :], modes[order, :], copy(nw), copy(cw),
        numeric_scale, resolved_gamma, objective, iterations, converged, Int(nstart), copy(excluded))
end

function _functional_fit(table, selected; k, species_col, trait_cols, categorical_traits,
    missing_traits, numeric_scale, trait_weights, gamma, nstart, maxiter, rng,
)
    for (value, name) in ((k, "k"), (nstart, "nstart"), (maxiter, "maxiter"))
        value isa Integer && !(value isa Bool) && 1 <= value <= typemax(Int) ||
            throw(ArgumentError("$name must be a positive Int-sized integer"))
    end
    work = _functional_workspace(table, selected; species_col, trait_cols, categorical_traits,
        missing_traits, numeric_scale, trait_weights, gamma)
    return _functional_fit_prepared(work; k, nstart, maxiter, rng)
end

"""
    functional_groups(traits::AbstractDataFrame; k, species_col=:species,
        trait_cols=nothing, categorical_traits=(), missing_traits=:error,
        numeric_scale=:sd, trait_weights=nothing, gamma=:auto,
        nstart=10, maxiter=100, rng=Random.default_rng())
    functional_groups(data::CommunityData, traits; same keywords...)

Cluster species by k-prototypes. Numeric contributions are weighted squared
differences; categorical contributions are weighted mismatches multiplied by
`gamma`. Numeric traits are divided by sample SD (`:sd`, default), range, or
one (`:none`); constant traits use scale one. A first-value offset preserves
small differences around large shared values. `trait_weights` is a dictionary
of positive per-column weights. `gamma=:auto` balances the mean pairwise numeric
and categorical contributions; it is one if either contribution is zero.

Numeric columns are inferred from non-Bool real values; strings, symbols and
Bool are categorical. Use `categorical_traits` for numeric category codes.
Within a categorical column, numeric codes must have one concrete type; string
and symbol values share their exact text. All numeric computation uses Float64.
Missing rows are errors or explicitly skipped. Data input selects its entire
species column domain, including unobserved species; table input fits all rows.
Extra table species are ignored only in the data overload. No input is mutated.

Starts use distinct randomly selected species and the supplied RNG. Empty
groups receive the farthest species from a group with more than one member.
Assignments retain their previous group on exact ties, and mode ties choose the
lexicographically smallest value. Results are numbered by first species. The
lowest objective among `nstart` runs is returned; `converged` and `iterations`
report the selected run. This local optimizer does not guarantee a global
minimum. Supply e.g. `StableRNG(42)` for reproducible random starts.

`numeric_centers` are in coordinates `(trait - numeric_offsets) / numeric_scales`;
result arrays are independently owned and should be treated as read-only.
"""
function functional_groups(traits::AbstractDataFrame; k, species_col = :species,
    trait_cols = nothing, categorical_traits = (), missing_traits::Symbol = :error,
    numeric_scale::Symbol = :sd, trait_weights = nothing, gamma = :auto,
    nstart = 10, maxiter = 100, rng::AbstractRNG = Random.default_rng(),
)
    return _functional_fit(traits, nothing; k, species_col, trait_cols, categorical_traits,
        missing_traits, numeric_scale, trait_weights, gamma, nstart, maxiter, rng)
end

function functional_groups(data::CommunityData, traits::AbstractDataFrame; k, species_col = :species,
    trait_cols = nothing, categorical_traits = (), missing_traits::Symbol = :error,
    numeric_scale::Symbol = :sd, trait_weights = nothing, gamma = :auto,
    nstart = 10, maxiter = 100, rng::AbstractRNG = Random.default_rng(),
)
    _validate_community_data(data)
    return _functional_fit(traits, data.species_labels; k, species_col, trait_cols, categorical_traits,
        missing_traits, numeric_scale, trait_weights, gamma, nstart, maxiter, rng)
end

function _functional_mapping(grouping::FunctionalGroupingResult)
    labels = _validated_labels(grouping.species_labels, length(grouping.membership), "species")
    groups = _validated_labels(grouping.group_labels, length(grouping.group_labels), "group")
    !isempty(groups) && length(grouping.group_sizes) == length(groups) ||
        throw(ArgumentError("functional grouping metadata dimensions do not match"))
    sizes = zeros(Int, length(groups))
    mapping = Dict{String, String}()
    for (label, group) in zip(labels, grouping.membership)
        1 <= group <= length(groups) || throw(ArgumentError("functional group membership is out of range"))
        sizes[group] += 1
        mapping[label] = groups[group]
    end
    all(>(0), sizes) && sizes == grouping.group_sizes ||
        throw(ArgumentError("functional group sizes do not match membership"))
    return mapping
end

function _functional_copy(grouping::FunctionalGroupingResult)
    return FunctionalGroupingResult(copy(grouping.species_labels), copy(grouping.membership),
        copy(grouping.group_labels), copy(grouping.group_sizes), copy(grouping.numeric_traits),
        copy(grouping.categorical_traits), copy(grouping.numeric_offsets), copy(grouping.numeric_scales),
        copy(grouping.numeric_centers), copy(grouping.categorical_modes), copy(grouping.numeric_weights),
        copy(grouping.categorical_weights), grouping.numeric_scale, grouping.gamma, grouping.objective,
        grouping.iterations, grouping.converged, grouping.nstart, copy(grouping.excluded_species))
end

"""
    collapse_traits(data::CommunityData, grouping::FunctionalGroupingResult;
                    statistic=:presence, normalize=false, missing_traits=:error)

Aggregate a community using a reusable species grouping. `:presence` returns
0/1 per group; `:richness` counts present species; `:abundance` sums original
abundances. `normalize=true` returns row proportions after the selected
statistic, with empty rows zero. The output is CSC `CommunityData` retaining
every fitted group. Unknown community species, including zero columns, are
errors unless `missing_traits=:skip` is explicit. Matching is exact.
"""
function collapse_traits(data::CommunityData, grouping::FunctionalGroupingResult;
    statistic::Symbol = :presence, normalize::Bool = false, missing_traits::Symbol = :error,
)
    _validate_community_data(data)
    _community_missing_mapping(missing_traits, "missing_traits")
    statistic in (:presence, :richness, :abundance) ||
        throw(ArgumentError("statistic must be :presence, :richness or :abundance"))
    mapping = _functional_mapping(grouping)
    unknown = filter(label -> !haskey(mapping, label), data.species_labels)
    isempty(unknown) || missing_traits === :skip ||
        throw(ArgumentError("community species missing from functional grouping: $(join(first(unknown, 5), ", "))"))
    result = aggregate_species(data, mapping; group_labels = grouping.group_labels,
        abundance_weighted = statistic === :abundance, missing_groups = :skip)
    statistic === :presence && fill!(nonzeros(result.matrix), 1)
    normalize || return result
    return CommunityData(_community_normalize_groups(result.matrix), result.sample_labels,
        result.species_labels, Val(:validated))
end

function collapse_traits(matrix::AbstractMatrix{<:Real}, grouping::FunctionalGroupingResult;
    species_labels = grouping.species_labels, sample_labels = nothing, kwargs...,
)
    data = CommunityData(matrix; species_labels, sample_labels)
    return collapse_traits(data, grouping; kwargs...)
end

"""
    functional_beta(data::CommunityData, grouping::FunctionalGroupingResult;
                    index_family=:sorensen, pairs=nothing, statistic=:presence,
                    missing_traits=:error)
    functional_beta(data::CommunityData, traits::AbstractDataFrame; k, ...)

Aggregate species into functional groups and apply [`beta_diss`](@ref). The
result contains the aggregated CSC `community`, an independent `grouping` and
the `beta::BetaDissResult`. `excluded_species` lists input species skipped for
this community; the grouping also records exclusions when it was fitted.
`statistic=:richness` or `:abundance` changes the
returned community; beta always uses presence of groups. Selected `pairs`,
empty-site NaNs and index families follow `beta_diss`. The table overload fits
`functional_groups(data, traits; k, ...)` first. A prepared grouping may be
reused for other communities without re-running clustering.
"""
function functional_beta(data::CommunityData, grouping::FunctionalGroupingResult;
    index_family::Symbol = :sorensen, pairs = nothing, statistic::Symbol = :presence,
    missing_traits::Symbol = :error,
)
    community = collapse_traits(data, grouping; statistic, missing_traits)
    beta = beta_diss(community; index_family, pairs)
    mapping = _functional_mapping(grouping)
    excluded = [label for label in data.species_labels if !haskey(mapping, label)]
    return FunctionalBetaResult(community, _functional_copy(grouping), beta, statistic, excluded)
end

function functional_beta(data::CommunityData, traits::AbstractDataFrame; k,
    species_col = :species, trait_cols = nothing, categorical_traits = (),
    missing_traits::Symbol = :error, numeric_scale::Symbol = :sd,
    trait_weights = nothing, gamma = :auto, nstart = 10, maxiter = 100,
    rng::AbstractRNG = Random.default_rng(), index_family::Symbol = :sorensen,
    pairs = nothing, statistic::Symbol = :presence,
)
    grouping = functional_groups(data, traits; k, species_col, trait_cols, categorical_traits,
        missing_traits, numeric_scale, trait_weights, gamma, nstart, maxiter, rng)
    return functional_beta(data, grouping; index_family, pairs, statistic, missing_traits)
end

function functional_beta(matrix::AbstractMatrix{<:Real}, input; species_labels = nothing,
    sample_labels = nothing, kwargs...,
)
    labels = species_labels === nothing && input isa FunctionalGroupingResult ? input.species_labels : species_labels
    data = CommunityData(matrix; species_labels = labels, sample_labels)
    return functional_beta(data, input; kwargs...)
end
