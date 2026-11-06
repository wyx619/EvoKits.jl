# Traits and community aggregation

[Documentation](README.md) · [Home](../README.md) · [简体中文](aggregation_zh.md)

On this page: [Trait summaries](#community-traits-and-species-groups) · [Functional groups](#mixed-traits-and-functional-beta) · [Group-count selection](#choose-a-functional-group-count) · [Lineages](#tree-groups-and-lineage-communities)

Summarize community traits and aggregate species into predefined, functional or lineage groups.

## Community traits and species groups

`map_trait` summarizes one numeric trait per site, and `aggregate_species`
combines species into predefined groups. Both accept a dense/CSC matrix or
`CommunityData`. Supply a vector in species-column order, a species-keyed
dictionary, or a `DataFrame`; table species must be unique and names match
exactly. Matrix calls can set `species_labels` and `sample_labels`.

```julia
using EvoKits, SparseArrays, DataFrames

data = CommunityData(sparse([2 1 0; 0 3 1; 0 0 0]);
    species_labels = ["A", "B", "C"], sample_labels = ["AB", "BC", "empty"])
traits = DataFrame(species = ["C", "A", "B"], trait = [-2.0, 2.0, 6.0])
summary = map_trait(data, traits; statistic = :mean, abundance_weighted = true)
summary.values # [10/3, 4, NaN]
summary.matched_richness # [2, 2, 0]

groups = Dict("A" => "group_1", "B" => "group_1", "C" => "group_2")
grouped = aggregate_species(data, groups) # species counts per group
proportions = aggregate_species(data, groups;
    abundance_weighted = true, normalize = true)
beta_diss(grouped) # grouped CommunityData works with taxonomic metrics
```

`map_trait` returns `CommunityTraitResult`: `values`, `sample_labels`, input
`richness`, used `matched_richness`, used `weight_sum`, and option metadata.
`statistic` is `:sum` (default), `:mean`, or `:sd`. By default every present
species contributes once; `abundance_weighted=true` uses original abundances.
Weighted SD uses reliability weights: the corrected denominator is
`sum(w) - sum(w^2)/sum(w)`, reducing to sample SD for unit weights. Set
`corrected=false` for population SD. Weighted mean/SD are invariant to common
weight scaling; integer abundances are not expanded into repeated records.
Numeric summaries use Float64; overflowing arithmetic is rejected. Summation
compensates for cancellation, and SD centers traits before online updates.

Missing traits or mappings for present species fail by default. Set
`missing_traits=:skip` or `missing_groups=:skip` to omit them. Unobserved
species may be unmapped; extra mapping species are ignored. Missing numeric
traits use Julia `missing`; matched NaN/Inf or nonnumeric traits are errors.
Empty/unmatched sites retain their labels and have zero sum, NaN mean/SD,
and zero used weight. Corrected SD with one matched species is NaN;
population SD is zero.

`aggregate_species` returns a new `CommunityData` whose CSC columns are groups.
Use [`community_longtable`](community_input.md#export-a-long-table-and-reconstruct-the-community)
to export grouped abundances or proportions with their complete label domains;
the same export accepts `collapse_traits` output and `collapse_range(...).community`.
By default entries count present species; `abundance_weighted=true` sums
abundance and preserves numeric type (Bool becomes Int). `normalize=true`
produces Float64 row proportions over retained groups; empty rows stay zero.
Group names are strings or integer identifiers. Columns follow first occurrence
in community species order, including mapped zero columns; `group_labels`
sets a complete order and can retain extra empty groups. Table selectors are
`species`/`trait` for `map_trait` and `species`/`group` for `aggregate_species`.
Inputs are not mutated, and CSC input is never expanded into a dense community
or a long occurrence table. Groups must be supplied; clustering and geographic
mapping are separate operations.

## Mixed traits and functional beta

`functional_groups` clusters numeric and categorical species traits using
k-prototypes. `collapse_traits` aggregates a community by that reusable mapping;
`functional_beta` applies the existing Sørensen or Jaccard decomposition to
presence of the resulting functional groups.

```julia
using EvoKits, DataFrames, SparseArrays, Random

data = CommunityData(sparse([1 1 0 0; 0 0 1 1; 1 0 1 0]);
    species_labels = ["A", "B", "C", "D"], sample_labels = ["west", "east", "mixed"])
traits = DataFrame(species = ["C", "A", "D", "B"],
    height = [10.0, 0.0, 12.0, 2.0], form = ["b", "a", "b", "a"])

groups = functional_groups(data, traits; k = 2, numeric_scale = :none,
    gamma = 10, nstart = 10, rng = MersenneTwister(42))
groups.membership            # [1, 1, 2, 2], in data species-column order
counts = collapse_traits(data, groups; statistic = :richness)
fb = functional_beta(data, groups)
fb.community.matrix          # sparse functional-group community
fb.beta.total                # pairwise functional beta dissimilarity
regions = phyloregion(fb.beta; k = 2)

# Fit and calculate in one call.
fb = functional_beta(data, traits; k = 2, rng = MersenneTwister(42))
```

The objective sums weighted squared numeric differences and weighted categorical
mismatches multiplied by `gamma`. Numeric traits use sample SD by default;
`numeric_scale=:range` or `:none` are also available, with constant traits
given scale one. Set positive per-column weights with e.g.
`trait_weights=Dict(:height => 2.0, :form => 0.5)`. The default `gamma=:auto`
balances mean numeric and categorical contributions over all species pairs;
it is one when either contribution is zero. Strings, symbols, Bool, and numeric
codes selected by `categorical_traits` are categorical. Numeric codes in one
categorical column must have the same concrete type.

Fitting a table uses all its species. Fitting with `data` selects the complete
community species-column domain, including all-zero columns. Names match
exactly; missing mappings or traits fail unless `missing_traits=:skip` is
explicit, in which case results report `excluded_species`. Select traits with
`trait_cols`. Numeric centers use `(value - numeric_offsets) / numeric_scales`.

Each start samples distinct species. The lowest-objective start is returned;
empty groups take a farthest species from a group with multiple members, and
ties in categorical modes use lexicographic order. Check `converged` and
`iterations`: this is a local optimizer and need not find the global minimum.
Pass a reproducibly seeded `rng` to reproduce starts and results.

`collapse_traits` supports `:presence`, `:richness`, `:abundance` and row
proportions, retaining every functional group and empty site. `functional_beta`
supports selected `pairs`; beta always uses group presence, even when the
returned community contains abundances. Scaling and initialization follow the
rules above and do not promise identical random clusters to R.

When comparing with R clustMixType, match scales, weights and initial
prototypes. Default weight estimation, ties and duplicate/empty-group handling
follow the conventions above.

### Choose a functional group count

`optimal_functional_groups` fits a collection of candidate counts and applies
the `lbow` heuristic underlying phyloregion's `quick_elbow`. It prepares species
matching, numeric scales, weights and gamma once so objectives are comparable.

```julia
selection = optimal_functional_groups(data, traits; k = 2:4,
    nstart = 20, rng = MersenneTwister(42))
selection.k
selection.objectives
selection.converged
selection.reason
if selection.grouping !== nothing
    selected_beta = functional_beta(data, selection.grouping)
end
```

Candidates must contain at least two distinct positive integers, no larger
than the complete-trait species count. They are fitted in ascending order,
consuming the supplied RNG sequentially across candidates and starts. Other
fitting keywords follow `functional_groups`; the table-only overload is also
available. The returned `FunctionalGroupSelectionResult` contains all candidate
objectives, iterations and convergence flags, plus `optimal_k` and the reusable
selected `grouping`. Only that grouping is retained in the result.

The normalized objective curve is evaluated with `low=0.08` and
`max_fraction=0.9`. This is a heuristic suggestion, without a claim of an
ecologically optimal count. Both recommendation fields are `nothing` if
`reason` is `:zero_objective`, `:nonconverged_candidates`,
`:nonmonotone_objective` or `:thresholds_not_met`; the diagnostics remain
available. Inspect the curve and convergence before using a recommendation.

## Tree groups and lineage communities

`get_clades` partitions a complete dated, ultrametric tree at an age before
the present, returning a `CladeGroupingResult` with tip-to-group membership,
group sizes, source CompactTree node IDs and cut metadata. Groups are numbered
by their first input tip; `tip_labels` and `membership` follow tree tip order.
No tree is pruned or changed.

```julia
mktempdir() do dir
    path = joinpath(dir, "tree.nwk")
    write_tree(path, "(((A:1,B:1):1,C:2):2,(D:2,E:2):2);")
    tree = read_tree(path)
    groups = get_clades(tree; cut = 2)
    groups.membership     # [1,1,1,2,2]
    groups.group_labels   # ["clade_1", "clade_2"]
    groups.group_nodes    # source node IDs, not tip positions

    data = CommunityData(sparse([2 4 0 0 0; 0 0 5 1 2; 0 0 0 0 0]);
        species_labels = tree.tip_labels, sample_labels = ["ab", "cde", "empty"])
    collapsed = collapse_range(data, groups) # reuse the full-tree grouping
    collapsed.community.matrix              # [1 0; 1 1; 0 0]
    counts = collapse_range(data, tree; cut = 2, statistic = :richness)
    counts.community.matrix                 # [2 0; 1 2; 0 0]
    mass = collapse_range(data, groups; statistic = :abundance)
    mass.community.matrix                   # [6 0; 5 3; 0 0]
    beta_diss(collapsed.community)           # direct downstream input
end
```

Supply exactly one of `cut` or target `k`. At a node exactly on an explicit
cut, `boundary=:older` (default) keeps the whole clade; `:younger` splits it.
At tree height, older gives one group. At cut=0, older can retain zero-age
clades; younger always gives singleton tips. Comparisons use represented ages
exactly, without an epsilon. For target k, simultaneous events and
multifurcations are split together; `actual_k` may exceed `requested_k`:

```julia
# For the tree above:
# groups = get_clades(tree; k = 3)
# groups.actual_k == 4, groups.cut == 2, groups.boundary == :younger
# get_clades(tree; cut=groups.cut, boundary=groups.boundary) reproduces the cut.
```

With k, boundary is chosen by the event; do not also specify `boundary`.
`k=1` keeps the root, and valid targets are 1:ntips. Unary nodes do not increase
the count. Ages use longest descendant paths, preserving short young branches
below long stems. Ultrametricity compares shortest/longest root-to-tip paths
with `atol=0.0`, `rtol=1e-8`; these tolerances validate the tree, not event ties.
Accepted near-ultrametric trees use longest descendant paths for node ages.
Branch lengths must be finite/nonnegative and path sums finite. Non-ultrametric
trees fail; date or transform them explicitly before calling. An explicit cut
uses O(nodes+tips) work/storage; target k also sorts branching events. No full
matrix of pairwise tip distances or descendant sets for every node is built.

`collapse_range` accepts a tree plus those cut keywords, or a prepared
`CladeGroupingResult`. It returns a `LineageCommunityResult` with labeled CSC
`community`, an independent `grouping`, `statistic` and `normalize` fields.
Default `:presence` is 0/1 per lineage; `:richness` counts present species;
`:abundance` sums original values and preserves their type (Bool becomes Int).
`normalize=true` computes row proportions after the chosen statistic, retaining
empty rows as zero; sums are checked for overflow. Matrix inputs use
`tip_labels` and `sample_labels`, defaulting to the complete tree/group tip order.

Name matching is exact. Missing tree species contribute zero, all full-tree
group columns are retained, and unknown species fail even if their column is
zero. Neither sparse input nor absent tree tips are expanded into dense padded
data, a dense membership matrix or a long occurrence table. Returned arrays
are independent; treat a reusable grouping's arrays as read-only. Group columns
are lineage identifiers: pass the collapsed community to `beta_diss`; the
original species tree cannot be used as a tree of these groups. This operation
does not truncate branch lengths as a tree time-slicing operation would.
