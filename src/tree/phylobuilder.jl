"""An independent completed tree and per-request taxonomic placement provenance."""
struct PhyloBuilderResult
    tree::CompactTree
    species_labels::Vector{String}
    status::Vector{Symbol}
    anchor_nodes::Vector{Union{Nothing, Int}}
    added_species::Vector{String}
    unmatched_species::Vector{String}
    extract::Bool
    branch_length::Union{Symbol, Float64}
end

function _phylobuilder_genus(label)
    parts = split(label, r"[ _]"; limit=2)
    return length(parts) == 2 && !isempty(parts[1]) && !isempty(parts[2]) ? String(parts[1]) : nothing
end

function _phylobuilder_records(species::AbstractVector; species_col, genus_col, taxonomy_cols)
    labels = _community_data_labels(species, length(species), "species")
    genera = [_phylobuilder_genus(label) for label in labels]
    return labels, genera, Pair{Symbol, Vector{Union{Nothing, String}}}[]
end

function _phylobuilder_records(table::AbstractDataFrame; species_col, genus_col, taxonomy_cols)
    String(species_col) in names(table) || throw(ArgumentError("taxonomy table lacks species column"))
    labels = _community_data_labels(table[!, species_col], nrow(table), "species")
    text(value) = ismissing(value) ? nothing : _community_identifier(value, "taxonomy")
    genera = String(genus_col) in names(table) ? text.(table[!, genus_col]) :
        [_phylobuilder_genus(label) for label in labels]
    (taxonomy_cols isa Tuple || taxonomy_cols isa AbstractVector) &&
        all(x -> x isa Union{Symbol, AbstractString}, taxonomy_cols) ||
        throw(ArgumentError("taxonomy_cols must list column names"))
    columns = String.(collect(taxonomy_cols))
    allunique(columns) && String(species_col) ∉ columns && String(genus_col) ∉ columns ||
        throw(ArgumentError("taxonomy_cols must be distinct from species and genus columns"))
    levels = Pair{Symbol, Vector{Union{Nothing, String}}}[]
    for column in columns
        column in names(table) || continue
        push!(levels, Symbol(column) => Union{Nothing, String}[text(value) for value in table[!, column]])
    end
    return labels, genera, levels
end

"""
    phylobuilder(species, tree; extract=true, species_col=:species, genus_col=:genus,
        taxonomy_cols=(:family,:suborder,:order,:class,:phylum),
        nonmonophyly=:error, branch_length=:ultrametric)

Explicitly complete a tree from exact species labels or a taxonomy table. The
table's species column contains full labels; genus and higher taxa are optional.
For vectors, infer genus only from labels containing an underscore or space.
Existing species stay in the tree. Missing species attach to the original
genus MRCA, or the parent when only one original congener exists. Named internal
nodes supply a fallback, first genus then the supplied taxonomy levels.

All anchors come from the original tree, so requests cannot seed later matches.
Nonmonophyletic genera fail by default; `nonmonophyly=:mrca` explicitly permits
their MRCA. Duplicate matching internal labels are ambiguous and fail. Unknown
taxa stay in `unmatched_species` and have status `:missing`. Duplicate requested
labels fail. This explicit operation does not change `align_community` defaults.

Default new terminal lengths are the anchor's longest descendant path, requiring
an ultrametric source (atol=0, rtol=1e-8). On other trees supply a finite
nonnegative terminal `branch_length`, or `:zero`. Original edges retain lengths.
`extract=true` then prunes to matched requested species; pruning may contract
unary paths and remove their labels. Node IDs/caches are rebuilt and all arrays
are independent. `anchor_nodes` always refers to original CompactTree IDs.
"""
function phylobuilder(species::Union{AbstractVector, AbstractDataFrame}, tree::CompactTree;
    extract::Bool=true, species_col=:species, genus_col=:genus,
    taxonomy_cols=(:family, :suborder, :order, :class, :phylum),
    nonmonophyly::Symbol=:error, branch_length=:ultrametric,
)
    nonmonophyly in (:error, :mrca) || throw(ArgumentError("nonmonophyly must be :error or :mrca"))
    mode = branch_length isa Symbol ? branch_length : _clade_time_value(branch_length, "branch_length")
    mode isa Symbol && mode ∉ (:ultrametric, :zero) && throw(ArgumentError("invalid branch_length mode"))
    _community_data_labels(tree.tip_labels, tree.ntips, "species")
    tree.ntips > 0 || throw(ArgumentError("source tree must contain tips"))
    all(x -> isfinite(x) && x >= 0, tree.edge_length) || throw(ArgumentError("invalid source branch lengths"))
    labels, genera, levels = _phylobuilder_records(species; species_col, genus_col, taxonomy_cols)
    isempty(labels) && throw(ArgumentError("species list must not be empty"))
    records = Dict(label => i for (i, label) in enumerate(labels))
    genus_tips = Dict{String, Vector{Int32}}()
    for (label, node) in zip(tree.tip_labels, tree.tip_ids)
        record = get(records, label, 0)
        genus = record == 0 ? _phylobuilder_genus(label) : genera[record]
        genus === nothing && continue
        push!(get!(genus_tips, genus, Int32[]), node)
    end
    named = Dict{String, Vector{Int32}}()
    for node in tree.postorder_internal
        label = tree.node_labels[node]
        isempty(label) || push!(get!(named, label, Int32[]), node)
    end
    index = _build_tree_query_index(tree)
    descendants = _tools_descendant_tip_counts(tree)
    genus_anchors = Dict{String, Int32}()
    statuses = fill(:missing, length(labels))
    anchors = Vector{Union{Nothing, Int}}(nothing, length(labels))
    for i in eachindex(labels)
        label, genus = labels[i], genera[i]
        if haskey(tree.tipname_to_id, label)
            statuses[i], anchors[i] = :tree, Int(tree.tipname_to_id[label])
            continue
        end
        if genus !== nothing && haskey(genus_tips, genus)
            tips = genus_tips[genus]
            anchor = get!(genus_anchors, genus) do
                if length(tips) == 1
                    parent = tree.parent_of_node[only(tips)]
                    return parent == 0 ? tree.root : parent
                end
                ancestor = foldl((a, b) -> _tree_lca(index, a, b), tips)
                nonmonophyly === :error && descendants[ancestor] != length(tips) &&
                    throw(ArgumentError("genus $genus is nonmonophyletic; choose nonmonophyly=:mrca explicitly"))
                ancestor
            end
            statuses[i], anchors[i] = :genus, Int(anchor)
            continue
        end
        choices = Pair{Symbol, Union{Nothing, String}}[:genus => genus]
        append!(choices, [level => values[i] for (level, values) in levels])
        for (level, taxon) in choices
            taxon === nothing && continue
            nodes = get(named, taxon, Int32[])
            isempty(nodes) && continue
            length(nodes) == 1 || throw(ArgumentError("ambiguous internal taxon label: $taxon"))
            statuses[i], anchors[i] = level, Int(only(nodes))
            break
        end
    end
    added = findall(i -> statuses[i] ∉ (:tree, :missing), eachindex(labels))
    ages = !isempty(added) && mode === :ultrametric ? _clade_node_ages(tree, 0.0, 1e-8) : nothing
    parents = copy(tree.parent_of_node)
    children = [copy(nodes) for nodes in tree.children]
    lengths = _edge_lengths_to_node(tree)
    tips = copy(tree.is_tip)
    node_labels = copy(tree.node_labels)
    tip_labels = copy(tree.tip_labels)
    single_root = nothing
    for i in added
        anchor = anchors[i]::Int
        if tips[anchor]
            # A single-tip source has no parent; place a zero-length root above it.
            if single_root === nothing
                single_root = Int32(length(parents) + 1)
                parents[anchor] = single_root
                push!(parents, 0)
                push!(children, Int32[anchor])
                push!(lengths, 0.0)
                push!(tips, false)
                push!(node_labels, "")
            end
            anchor = Int(single_root)
        end
        terminal = mode === :ultrametric ? ages[anchors[i]::Int] : mode === :zero ? 0.0 : mode
        node = Int32(length(parents) + 1)
        push!(parents, Int32(anchor))
        push!(children, Int32[])
        push!(children[anchor], node)
        push!(lengths, terminal)
        push!(tips, true)
        push!(node_labels, labels[i])
        push!(tip_labels, labels[i])
    end
    # The low-level constructor associates labels with ascending tip node IDs.
    built = _compact_tree_from_topology(parents, children, lengths, tips, tip_labels, node_labels)
    if extract
        matched = labels[statuses .!= :missing]
        isempty(matched) && throw(ArgumentError("no requested species can be matched or placed"))
        built = keep_tip(built, matched)
    end
    all(isfinite, built.dist_from_root) || throw(ArgumentError("completed root-to-tip paths overflowed"))
    return PhyloBuilderResult(built, labels, statuses, anchors, labels[added],
        labels[statuses .== :missing], extract, mode)
end
