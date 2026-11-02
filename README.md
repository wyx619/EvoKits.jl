# EvoKits.jl

EvoKits is a Julia-native toolkit for practical phylogenetic work. It provides
tree I/O, a compact in-memory tree representation, topology utilities,
community phylogenetic metrics, tree identity mapping, and lightweight
comparative analysis.

## Install and Load

```julia
using Pkg
Pkg.add(url = "https://github.com/wyx619/EvoKits.jl")
using EvoKits
```

EvoKits loads MKL as its first dependency and supports selecting Intel MKL as
the BLAS backend.

## Tree I/O and Topology

```julia
using EvoKits

tree = read_tree("example.nwk")

ultrametric_tree = make_ultrametric(tree)
binary_tree = multi2di(ultrametric_tree; random = false)
write_tree("binary.nwk", binary_tree)
```

`CompactTree` is the common input type for EvoKits tools. It caches topology,
edge lengths, traversal order, and tip-label lookup. `read_tree` reads a Newick
file directly into this representation.

Tree queries use compact node IDs, tip positions, or tip labels explicitly:

```julia
mrca(tree, "A", "B")
distance(tree, 1, 2)
pairwise_distance(tree)
descendant_tips(tree, tree.root; labels = true)
tree_imbalance(tree; metric = :sackin)
```

Weighted distances use branch lengths; pass `weighted = false` for edge
counts. The source is organized into four layers: `base`, `tree`, `traits`,
and `community`.

## Prune and Map Trees

```julia
pruned = keep_tip(tree, ["A", "C", "D"])
without_b = drop_tip(tree, "B")

map = build_phyloref(pruned)
node_map = phyloref_node_table(pruned; map = map)
edge_map = phyloref_edge_table(pruned; map = map)

ape_edges = ape_edge_matrix(pruned; map = map)
signature = phylo_edge_signature(pruned, 1; map = map)
```

`PhyloRef` provides stable node and edge identities, including conversions
between `CompactTree` IDs and R/ape numbering. Its tables and signatures make
tree results easier to compare across tools and runtimes.

## Community Metrics

`CommunityData` keeps community values and row/species labels together. It
accepts a matrix, a long `DataFrame`, or indexed `(rows, cols, values)` vectors.
Matrix inputs are copied and retain dense/CSC storage and numeric type;
`storage=:sparse` requests CSC storage. Long tables and triplets build sparse
matrices directly, including when the input is too large for a dense table.

```julia
using DataFrames

records = DataFrame(
    sample = ["s1", "s1", "s2", "empty"],
    species = ["A", "C", "A", "A"],
    value = [2, 1, 3, 0],
)
data = CommunityData(records; sample=:sample, species=:species, value=:value)
data.matrix
data.sample_labels   # ["s1", "s2", "empty"]
data.species_labels  # ["A", "C"]

weighted_endemism(data).we
beta_diss(data; pairs=[(1, 2)])
phylo_endemism(data, tree).pe

# Align once to the complete tree, or also cache branch sets and occupancy.
aligned = align_community(data, tree)
ctx = prepare_community(aligned, tree)
pd(ctx).pd
```

Labels default to first-occurrence order for long tables. Explicit
`sample_labels` and `species_labels` set a complete output order and can retain
empty sites or unobserved species. String and integer identifiers are accepted;
integers become strings. Labels must be unique and nonempty. Matrix labels
default to `sample_1`, ... and `species_1`, ...; supply actual species names for
tree-based analysis. Record values must be finite, non-negative real numbers
without missing values. `value=nothing` counts occurrences instead of reading
an abundance column.

Duplicate records sum by default; `duplicates=:error` rejects them, including
zero-valued duplicates. Overflowing sums fail. Values retain their abundance
magnitudes; presence-based metrics use every positive value as presence.
For example, indexed triplets can be constructed as:

```julia
data = CommunityData([1, 1, 2], [1, 2, 1], [2, 1, 3];
    sample_labels=["s1", "s2"], species_labels=["A", "C"])
```

Tree-based calls align species by name, fill unrepresented tree species with
zero columns, and preserve the complete tree and its root. Unknown input
species always cause an error. `align_community(data, tree; missing_tips=:error)`
also rejects missing tree species. Alignment preserves dense/CSC storage;
aligning once lets subsequent calls reuse the aligned matrix. Prepared contexts
add branch and occupancy reuse. Treat `CommunityData` arrays as read-only;
constructors, explicit alignment, metric results, and contexts own their arrays.

All existing matrix-based community metrics and `prepare_community` accept
`CommunityData`. Results keep input row order; existing labeled result types
receive its site labels. Set labels at construction instead of supplying
`tip_labels` or `sample_labels` to metric calls. SES uses sparse species sets
for observations and null draws, including abundance-weighted MPD/MNTD. Several
other legacy metrics still require dense internal workspaces. Existing direct
matrix calls remain available:

```julia
community = [
    1 0 1
    0 1 1
]

pd_result = pd(community, tree; tip_labels = tree.tip_labels)
ses_result = ses_mpd(community, tree; runs = 999)
```

Available community tools include phylogenetic diversity (`pd`), evolutionary
distinctiveness (`evol_distinct`), and standardized-effect-size metrics for
PD, MPD, and MNTD (`ses_pd`, `ses_mpd`, and `ses_mntd`).

Community beta diversity and diversity extensions use the same row/column
contract:

```julia
comdist(community, tree)
comdistnt(community, tree)
unifrac(community, tree)
phylosor(community, tree)
raoD(community, tree)

psv_result = psv(community, tree)
psv_result.psv
expected_pd(tree).expected_pd
```

`raoD` returns a `RaoDResult` containing `Dkk`, `Dkl`, `H`, `alpha`, `beta`,
and `Fst`, matching the Picante partitioning convention. It requires an
ultrametric tree and uses distance-to-MRCA (`patristic / 2`).

Rows are communities and columns are tips. Pass `tip_labels` when input
columns are not already in tree-tip order. Empty communities return `NaN`
where a distance or diversity value is undefined; singleton communities are
accepted by PD, UniFrac, PhyloSor, and RaoD.

The non-abundance-weighted MPD/MNTD paths use a virtual-tree calculation, so a
community query does not require scanning every edge of a large source tree.
Null means and standard deviations are accumulated with Welford's online
algorithm.

`ses_pd`, `ses_mpd`, and `ses_mntd` keep only one null replicate and reusable
species/tree workspaces, with no dense site-by-species randomization buffers.
Supported null models are `:taxa_labels`, `:richness`, `:frequency`,
`:sample_pool`, `:phylogeny_pool`, `:trialswap`, and `:independentswap`.
The storage refactor retains their existing draw rules and permutation-based
RNG consumption: richness preserves row sums, frequency preserves column sums,
and swaps preserve both. The existing `:sample_pool` path uses the same full
tree species pool as `:richness`; `:phylogeny_pool` additionally shuffles labels.
Abundance-weighted MPD/MNTD randomize presence while keeping original abundances
at their tip positions. Zero-weight null tips remain MNTD neighbor candidates.

All three SES methods use finite draws for null means and sample SD. No finite
draws give `NaN` mean/SD/z; a single finite draw gives a mean and `NaN` SD/z;
zero null SD gives `NaN` z. Rank p-values retain the existing convention
`(1 + count_below + half_ties) / (runs + 1)`; nonfinite draws contribute no rank
increment. PD null draws now use a fresh branch-visit generation each time,
fixing skipped branches between repeats. Without the root, PD sums below the
MRCA directly and avoids subtracting a potentially long shared root path.
SES virtual edges also fall back to summing original branches when depth
subtraction would lose short distances beneath a long shared path.

For Caper-compatible ED handling on polytomies, use
`evol_distinct(tree; type = :fair_proportion,
polytomy_correction = :isaac)`. The other supported modes are `:mooers` and
`:none`.

## Endemism, Relative PD and Beta Decomposition

These functions accept dense matrices or `SparseMatrixCSC` inputs. Rows are
sites, columns are species, and every positive value denotes presence. Values
must be finite and non-negative. For direct matrix calls to tree-based functions, supply one column
per tree tip, including zero columns
for species absent from the collection. Use `tip_labels` when the column order
differs from `tree.tip_labels`. `CommunityData` calls perform this alignment
and zero padding automatically.

```julia
using SparseArrays

comm = sparse(community)
we = weighted_endemism(comm)
pe = phylo_endemism(comm, tree; tip_labels = tree.tip_labels)
strict_pe = phylo_endemism(comm, tree; weighted = false)
rpd = relative_pd(comm, tree)

we.we
pe.pe
strict_pe.pe
rpd.rpd
rpd.pd_normalized
rpd.pd_equal_normalized

# Taxonomic beta uses species counts and does not require a tree.
tax_beta = beta_diss(comm; index_family = :sorensen)
tax_beta.total
tax_beta.turnover
tax_beta.nestedness
tax_selected = beta_diss(comm; index_family = :jaccard, pairs = [(1, 2)])

beta = phylobeta(comm, tree; index_family = :sorensen)
beta.total
beta.turnover
beta.nestedness

# Calculate selected pairs without allocating three full site-by-site matrices.
selected = phylobeta(comm, tree; index_family = :jaccard, pairs = [(1, 2)])
selected.total
```

Weighted endemism (WE) sums reciprocal species occupancy. Weighted phylogenetic
endemism (PE) sums branch length divided by the number of sites containing that
branch; strict PE sums lengths of branches found in exactly one site. Occupancy
uses the complete input collection, so independently calculating PE or WE on
separate batches changes the statistic. Empty sites return zero.

Relative PD (RPD) compares the original tree with the same topology after
equalizing positive branch lengths and retaining zero branches. Each tree is
normalized to total length one before calculating the ratio. All tree-based
metrics in this section include paths to the supplied root; the input tree is
preserved. RPD returns `NaN` for zero denominators.

`beta_diss` returns a `BetaDissResult` partitioning taxonomic beta diversity
using shared and unique species counts. It requires no tree and ignores
abundance magnitudes. `phylobeta` uses shared and unique branch lengths.
Both support Sørensen and Jaccard decomposition. Their `total`, `turnover`
and `nestedness` fields are symmetric matrices, or vectors in the supplied pair
order when `pairs` is used. For defined results, total equals turnover plus
nestedness. Pair requests retain repeated, reversed, and self pairs.

For `beta_diss`, comparing an empty and a nonempty site gives total one and
`NaN` turnover/nestedness. Two empty sites give `NaN` for all components,
including the diagonal; nonempty diagonal entries are zero. For `phylobeta`,
comparing a zero-PD site with a positive-PD site gives total one
and undefined (`NaN`) turnover/nestedness; two zero-PD sites give `NaN` for all
three components. Positive-PD diagonal entries are zero.

`weighted_endemism`, `phylo_endemism`, `relative_pd`, `beta_diss`, and `phylobeta` accept
`sample_labels` and retain row order. Results include richness or occupancy
diagnostics where applicable. `pd`, `unifrac`, and `phylosor` also use sparse
presence processing without constructing dense community inputs. WE, PE and PD
are computed without retaining all tip-to-root paths. Branch lengths must be
finite and non-negative, with a finite total branch length.

For repeated analyses of the same collection, prepare it once:

```julia
ctx = prepare_community(comm, tree; tip_labels = tree.tip_labels)
pd(ctx).pd
pd(ctx; include_root = false).pd
weighted_endemism(ctx).we
phylo_endemism(ctx).pe
phylo_endemism(ctx; weighted = false).pe
relative_pd(ctx).rpd
unifrac(ctx)
phylosor(ctx)

# Each request reuses the label alignment and branch sets.
first_batch = phylobeta(ctx; pairs = [(1, 2)])
second_batch = phylobeta(ctx; index_family = :jaccard, pairs = [(2, 1)])
```

`prepare_community` returns a `CommunityContext` with independent snapshots of
label alignment, per-site branch sets, global species/branch occupancy, and
per-site summaries. Pair indices always refer to the original input rows;
WE/PE occupancy always refers to the complete collection used to prepare the
context. Species occupancy retains the input column order. `sample_labels`
are set during preparation, and `pd_dataframe(ctx)` is also supported.

Changing the original matrix, tree, or labels does not update the context;
prepare a new one when the data changes. Returned results may be modified
independently. Treat the context's internal arrays as read-only. Preparation
retains branch sets and summaries, so it costs time and memory up front and is
intended for reuse. It does not retain the input tree or matrix, or cache full
pairwise matrices. Full matrix calls still perform pairwise calculations and
allocate their outputs; selected-pair calls allocate only requested results.

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

## Regionalization from beta diversity

`phyloregion` clusters a complete distance matrix and cuts the full hierarchy
into an explicit number of regions `k`, using UPGMA by default. It also
accepts full `BetaDissResult` and `PhylobetaResult` objects, retaining site labels.
Select `component=:total` (default), `:turnover`, or `:nestedness`.

```julia
regional_data = CommunityData(sparse([1 1 0 0; 1 1 1 0; 0 0 1 1; 0 1 1 1]);
    species_labels = ["A", "B", "C", "D"],
    sample_labels = ["west_1", "west_2", "east_1", "east_2"])
regions = phyloregion(beta_diss(regional_data); k = 2)
regions.membership       # input-row order; region IDs start at 1
regions.cluster_sizes
regions.region_distances
regions.distinctiveness

# Or supply a complete distance matrix and its labels.
regions = phyloregion([0 1 8 10; 1 0 6 12; 8 6 0 2; 10 12 2 0];
    k = 2, sample_labels = ["A", "B", "C", "D"])
# membership = [1, 1, 2, 2]; region_distances = [0 9; 9 0]
```

`PhyloregionResult` returns site `sample_labels` and `membership`,
`cluster_sizes`, `region_distances`, `distinctiveness`, and the complete
hierarchy (`merge` and `heights`), along with `k`, `method`, and source
`component`. Regions are numbered by their first input site. Negative merge
entries identify input rows; positive entries identify earlier merge rows.
Region distances average all cross-site pairs between two regions and have
zero diagonals. Distinctiveness averages distances to the other regions with
equal region weight, regardless of region size; it is NaN for k=1.

Input must be nonempty, square, exactly symmetric, finite, non-negative, with
zero diagonal; it need not be Euclidean. Float64 conversion must preserve
finite positive distances. Undefined values, including NaN diagonals from
empty sites, are rejected; explicitly select a defined subset and its labels
before calling. No sites are automatically removed or distances imputed.
Selected-pair beta results cannot supply the complete distances needed here.

Supported methods follow R `hclust` recurrences:

| `method` | Meaning |
|---|---|
| `:average` | UPGMA; average of all cross-site distances |
| `:single` | Minimum cross-site distance |
| `:complete` | Maximum cross-site distance |
| `:mcquitty` | WPGMA; equal weight for the two child clusters |
| `:ward_D` | R's older Ward recurrence on supplied distances |
| `:ward_D2` | Ward recurrence on squared distances, square-root heights |
| `:centroid` | UPGMC recurrence on supplied distances |
| `:median` | WPGMC recurrence on supplied distances |

`Symbol("ward.D")` and `Symbol("ward.D2")` are accepted aliases; results use
`:ward_D` and `:ward_D2`. Ward's variance interpretation requires Euclidean
distances. Centroid/median do not automatically square the input; their heights
can decrease and remain in chronological merge order. A linkage update outside
finite Float64 range is an error; rescale input distances in that case.

The six monotone methods use O(n²) nearest-neighbor chaining. Centroid/median
use global minima with cached row minima, with O(n³) worst-case time. Both use
a packed Float64 work triangle of `n*(n-1)/2` entries. The full input and
regional output still require quadratic memory. Ties are deterministic for
fixed row order: chains start at the earliest active cluster, retain the
predecessor on equal nearest distances, otherwise choose the earliest active
cluster; equal-height discoveries keep their order. Centroid/median choose the
lexicographically first active pair. Ties can differ from R `hclust` and depend
on row order. Inputs are not mutated; outputs own their arrays.

### Reuse a hierarchy and evaluate cuts

```julia
ctx = prepare_phyloregion(beta_diss(regional_data); method = :average)
cutree(ctx; k = 2)                 # memberships only; O(n), no reclustering
regions3 = phyloregion(ctx; k = 3) # regional summaries, reusing the same tree
cutree(regions3; k = 2)            # any PhyloregionResult retains the full tree

selection = select_linkage(ctx)    # compare all eight methods
selection.methods
selection.correlations            # original vs cophenetic Pearson correlation
selection.best_method             # first finite maximum, or nothing

evaluation = optimal_phyloregion(ctx; k = 4) # evaluate every cut from 1 to 4
evaluation.df                     # columns: k, ev, totwss, totbss, tss
evaluation.optimal_k               # suggested k, or nothing
evaluation.reason
```

`PhyloregionContext` owns an original packed distance snapshot and full tree;
input changes do not affect it. Treat context arrays as read-only. Preparation
retains one packed triangle and uses a second as temporary clustering
workspace. Repeated `cutree` calls only cut the tree; `phyloregion(ctx; k)` also
visits site pairs for regional means and allocates the k×k regional matrix.

`select_linkage` also accepts a distance matrix or full beta result. Restrict
methods with, for example, `methods=(:average, :complete, :single)`. It compares
original distances with the height where each pair first joins the tree.
Larger correlation means better distance preservation, not ecological
correctness. Constant distances/heights and a single pair give NaN; if all
scores are undefined, `best_method` is `nothing`. No full cophenetic matrix is
allocated. With a context, its existing tree is reused for the matching method.

`optimal_phyloregion` also accepts a distance matrix or full beta result,
building one hierarchy. `k` is the upper bound, defaulting to `min(n,20)`.
Within-cluster sums are the sum of squared unordered distances divided by the
cluster's site count; TSS uses all sites, BSS=TSS−WSS, and `ev=BSS/TSS`.
For non-Euclidean dissimilarities these are distance-based scores; EV may be
nonmonotone or outside [0,1].

Selection follows the saved phyloregion threshold heuristic: try
`ev_thresholds=(0.95,0.9,0.8,0.75,0.67,0.5,0.33)` descending, then
`inc_thresholds=(0.01,0.05,0.1)` ascending. Find the first k meeting each EV
threshold and the first whose next EV increment is small enough, and recommend
their maximum for the first successful threshold combination. Tolerance is
`10^(-precision)`, default precision=3; the last evaluated k has no next
increment. This is a heuristic, not a curvature optimization. Unmet thresholds
or zero dispersion return `optimal_k=nothing` and the full table, with
`:thresholds_not_met` or `:zero_dispersion` as the reason. Sums use original
squared-distance units and can overflow/underflow, while EV and decisions use
normalized distances. Scores visit all pairs once with O(n) extra storage.

These functions add no dependencies. NMDS, geographic constraints, spatial
polygons and plotting are outside this numerical API.

## NRI and NTI

For a single observed community and a local species pool, use `nri_nti`:

```julia
pool = ["A", "B", "C", "D"]
observed = ["A", "B", "D"]

result = nri_nti(tree, pool, observed; nreps = 999, seed = 20260713)
result.nri
result.nti
result.status   # :ok or :null_sd_zero
```

The pool is sampled without replacement at the observed richness. Both tip
labels and tip positions (`1:tree.ntips`) are accepted. `nri_nti` uses a stable
random generator and a local virtual-tree MPD/MNTD calculation; it does not
assume a global species pool or retain all null replicates in memory. The
reference field names `NRI`, `NTI`, `MPD_obs`, and `MNTD_obs` are also available
as compatibility properties.

For a batch of spatial grids, pass `grid_id` so the global seed is stably
derived per grid, independent of processing order. This is the convention used
by the validated Swq workflow:

```julia
result = nri_nti(tree, pool, observed;
    nreps = 999,
    seed = 20260713,
    grid_id = "1745",
)
```

The effective derived seed is stored in `result.seed`. Do not combine `grid_id`
with an explicit `rng`; an injected RNG already defines the random stream.

For many grids, `nri_nti_table` accepts two standard delimited long tables with
exactly the columns `grid_id,species`: one local-pool table and one
observed-community table. Tables may be CSV, TSV, TXT, or another delimiter
recognized by CSV.jl; use `delim = '\t'` when an explicit TSV delimiter is
needed. Each grid must be a contiguous block in each file. Grids are processed
in parallel when Julia starts with multiple threads.

```julia
batch = nri_nti_table(tree, "local_pool.tsv", "observed_community.tsv";
                        nreps = 999, seed = 20260713, threaded = true,
                        delim = '\t')
batch.grid_ids
batch.results
```

## Comparative Utility

### Pagel's lambda

`pagel_lambda` estimates Pagel's lambda for a continuous univariate trait under a
phylogenetic GLS model. The tree must be rooted, strictly bifurcating, and
ultrametric; trait values must be finite and supplied in tree-tip order.

```julia
trait = [1.2, 2.1, 1.8, 2.4]
lambda_result = pagel_lambda(tree, trait)
lambda_result.lambda
lambda_result.ci_lower, lambda_result.ci_upper
lambda_result.p_value
summary_pagel_lambda(lambda_result)
```

The estimate is constrained to `[0, 1]`: values near zero indicate little
phylogenetic covariance after accounting for the fitted residual variance,
whereas values near one retain the Brownian covariance structure. The result
is a `PagelLambdaResult` with `success`, `lambda`, `loglik`, `aic`, `sigma2`,
`n`, convergence counters, profile-likelihood interval fields
`ci_lower`/`ci_upper`, likelihood-ratio fields `lr_stat`/`p_value`, and the
underlying `fit` tuple. `summary_pagel_lambda` returns the same information
with conventional `logLik` and `AIC` field names.

For tabular data, pass the trait column and optionally a taxon column. The
taxon column is aligned to `tree.tip_labels` before fitting:

```julia
using DataFrames

data = DataFrame(species = tree.tip_labels, body_size = trait)
lambda_result = pagel_lambda(
    tree,
    data,
    :body_size;
    taxon_col = :species,
    profile_points = 101,
)
```

`root_model` accepts `:random` or `:fixed`; optimization controls are
`max_iterations`, `rel_tol`, and `profile_points`. Use `align_traits_to_tree`
when a vector or table must be reordered independently of model fitting.

### Phylogenetic K (PhyloK)

`phylogenetic_signal_k` implements the Picante/phytools K signal statistic for
continuous traits. It compares the observed PIC variance with a null
distribution generated by shuffling trait values across the tree tips:

```julia
using Random

k_result = phylogenetic_signal_k(
    tree,
    trait;
    permutations = 999,
    rng = MersenneTwister(20260717),
)

k_result.k
k_result.p_value
k_result.z_score
k_result.pic_variance_observed
```

The return value is a `PhyloSignalResult` containing `k`, observed/null PIC
variance summaries, `p_value`, `z_score`, reordered `trait` values, `tip_labels`,
`valid_replicates`, `permutations`, and `status`. `tip_labels` may be supplied
when the trait vector is not already in tree order:

```julia
k_result = phylogenetic_signal_k(
    tree,
    [1.8, 1.2, 2.1, 2.4];
    tip_labels = ["C", "A", "B", "D"],
    permutations = 999,
)
```

As a rule of thumb, `K approximately 1` is Brownian-like, `K > 1` indicates stronger
phylogenetic conservatism, and `K < 1` indicates weaker signal than expected
under the Brownian reference. The null p-value is a rank-based value from the
valid tip-shuffle replicates; pass an explicit `AbstractRNG` for reproducible
permutations. The default RNG is `StableRNG(20260713)`.

PhyloK currently requires finite traits and the same strictly bifurcating,
positive-branch-length tree conditions used by its PIC variance calculation.
Missing values are rejected. Set `permutations = 0` to obtain the observed
statistic without generating a null distribution.

The restricted P2 comparative utilities are independent contrasts and the
Picante/phytools K signal test:

```julia
pic = phylogenetic_contrasts(tree, trait)
signal = phylogenetic_signal_k(tree, trait; permutations = 999)
signal.k
signal.p_value
```

PIC currently requires strictly bifurcating trees with positive branch
lengths. Use `multi2di` explicitly when a multifurcating tree must be resolved;
missing trait values are rejected.

## Phylogenetic D

`phylo_d` measures phylogenetic signal in a binary trait using the
Fritz-Purvis D statistic. It compares the observed sum of nodal changes with
tip-state randomizations and Brownian-threshold simulations.

```julia
using Random

trait = [1, 1, 0, 0]  # tree.tip_labels order
d = phylo_d(tree, trait; permutations = 999, rng = MersenneTwister(42))

d.d             # D estimate; 0 is Brownian-like, 1 is random-like
d.p_random      # random-null probability
d.p_brownian    # Brownian-null probability
d.observed
d.mean_random
d.mean_brownian
```

The trait must contain exactly two finite states and no missing values. String
states and `Bool` values are accepted. For a DataFrame, use
`phylo_d(tree, data, :trait; taxon_col = :species)` to align rows by tip
labels. `random_weights` optionally supplies weighted sampling probabilities
for the random null model. Set `store_nodal_values = true` to retain the
permutation-by-node nodal value matrices; they are omitted by default to keep
large analyses memory-efficient. Compatibility properties `D`, `DEstimate`,
`Pval1`, `Pval0`, `Parameters`, `Permutations`, and `NodalVals` mirror the
original caper result names.

For clade-wise analyses, `phylo_d_subset` applies the same calculation to
filtered internal clades and returns one summary row per clade.

## EDGE2

EDGE2 combines evolutionary distinctiveness with species extinction
probabilities. It uses the same `CompactTree` representation as the other
EvoKits tools, so no separate tree parser or `Phylo.jl` dependency is needed.

```julia
pext = [0.10, 0.20, 0.30, 0.40]  # tree.tip_labels order
edge2 = edge2_mod(tree, pext)

edge2.species
edge2.edge       # EDGE2 scores
edge2.ed         # evolutionary distinctiveness
edge2.pd         # total phylogenetic diversity
edge2.epd_loss   # expected PD loss
```

The reference field spellings `EDGE`, `ED`, `TBL`, `PD`, and `ePDloss` remain
available as compatibility properties. Extinction probabilities must be
finite values in `(0, 1]` and must follow the tree tip order.

For category-based workflows, use `PextDistribution`,
`calc_pext_distribution`, `load_pext_distribution_csv`, and `sample_pext`.
`reorder_pext_to_tree` validates species coverage before returning values in
tree order.

## End-to-End Example

This example exercises the tree, community, and restricted trait layers in one
session:

```julia
using EvoKits
using Random

write_tree("example.nwk", "((A:1,B:1):1,(C:1,D:1):1);")
tree = read_tree("example.nwk")

community = [
    1 1 0 0
    0 1 1 0
]
pd_result = pd(community, tree)
beta = phylosor(community, tree)
rao = raoD(community, tree)

trait = [1.0, 2.0, 3.0, 4.0]
pic = phylogenetic_contrasts(tree, trait)
signal = phylogenetic_signal_k(tree, trait;
    permutations = 99, rng = MersenneTwister(7))

pd_result.pd
beta[1, 2]
rao.Dkk
pic.contrasts
signal.k
```

## Public API Contract

### Shared Ordering

`tree.tip_labels` is the canonical tip order. A tip position is an integer in
`1:tree.ntips`; a compact node ID is an integer in `1:tree.nnodes`. Functions
that accept labels resolve them through `tree.tipname_to_id`.

Community matrices use rows as communities and columns as tips. Row order is
preserved in every returned vector or matrix. `tip_labels` describes the input
column order and reorders columns into `tree.tip_labels` order; it does not
rename or reorder community rows. Duplicate, missing, or unknown labels raise
`ArgumentError`.

Community values must be finite and non-negative. Values greater than zero
indicate presence; abundance-weighted methods use positive values as weights.

Trait vectors are aligned to tree order by default. When `tip_labels` is
provided, it describes the order of the supplied trait vector and the vector
is reordered to `tree.tip_labels`. Missing and non-finite trait values are
rejected by K and PIC.

### Matrix Semantics

`comdist`, `comdistnt`, `unifrac`, and `phylosor` return ordinary `n x n`
`Float64` matrices. Matrix indices follow community row order.

- `comdist`: mean pairwise distance across the two communities. The diagonal
  is within-community mean distance; an empty community gives `NaN`.
- `comdistnt`: symmetric mean of the two directional nearest-taxon means. The
  diagonal is the within-community nearest-taxon mean. With
  `exclude_conspecifics=true`, a singleton diagonal is `NaN`.
- `unifrac`: unweighted branch-length distance. The diagonal is zero,
  including an empty-versus-empty comparison; an empty-versus-nonempty pair is
  one when the nonempty community spans positive branch length.
- `phylosor`: `2 * shared_PD / (PD1 + PD2)`. The diagonal is one for a
  nonempty community and zero for an empty community.

`raoD` returns `RaoDResult` with `Dkk`, `Dkl`, `H`, `total`, `alpha`, `beta`,
and `Fst`, matching the Picante partitioning convention. It requires an
ultrametric tree and uses distance-to-MRCA (`patristic / 2`).

### Community Result Types

- `pd` returns `PDResult` with `sample_labels`, `pd`, `richness`, and
  `include_root`.
- `weighted_endemism` returns `WeightedEndemismResult` with `sample_labels`,
  `we`, `richness`, and `species_occupancy` in input column order.
- `phylo_endemism` returns `PhyloEndemismResult` with `sample_labels`, `pe`,
  `richness`, `branch_occupancy` in tree edge order, and `weighted`.
- `relative_pd` returns `RelativePDResult` with `sample_labels`, `rpd`,
  `pd_normalized`, `pd_equal_normalized`, `richness`, `total_branch_length`,
  and `positive_branch_count`.
- `phylobeta` returns `PhylobetaResult` with `sample_labels`, `total`,
  `turnover`, `nestedness`, `index_family`, and `pairs`.
- `ses_pd`, `ses_mpd`, and `ses_mntd` return typed results containing observed
  values, null means/SDs, z scores, rank p-values, run count, and null model.
  Empty communities have `NaN` MPD/MNTD. Singleton observations are `NaN`,
  except abundance-weighted MPD with positive mass, which is zero.
- `psv` returns `PSVResult` with `psv`, `richness`, `variance`, and `scale_vcv`.
  Communities with fewer than two taxa have `NaN` PSV and variance. `psr`,
  `pse`, `psc`, and `psd` return named tuples whose field names are the metric
  names documented by their function.
- `expected_pd` and `pd_variance` return `ExpectedPDResult`; `richness` is
  `1:tree.ntips`, `expected_pd` is the analytical expectation, and
  `variance_pd` is the upper-bound or exact Picante variance selected by
  `upper_bound`.
- `pd_calc(tree, tips; variant=...)` returns one scalar for a tip subset.
  Matrix input returns a named tuple with `pd`, `richness`, and normalized
  `variant`. Accepted variants are `TBL`, `MST`, `UEH`, `SBL`, and `TIP`.
- `pd_bootstrap(tree, ntips; ...)` returns a vector of replicate scalars.
  `pd_bootstrap(community, tree; ...)` returns `PDBootstrapResult` with
  `observed`, `replicates`, and `richness`.

### Trait Result Types

`pagel_lambda` returns `PagelLambdaResult` with the fitted lambda, profile
interval, likelihood-ratio test, fit diagnostics, and model summary fields
described above. The DataFrame overload performs taxon-label alignment before
fitting.

`phylogenetic_signal_k` returns `PhyloSignalResult` with `k`, observed and
null PIC variance, null SD, rank p-value, z score, reordered trait values, tip
labels, valid replicate count, and status. Its null model is tip shuffling;
pass an explicit `rng` for reproducibility. The default is
`StableRNG(20260713)`.

`phylogenetic_contrasts` returns `PICResult` with `contrasts`, internal-node
IDs, expected contrast distances, local ancestral values, tip labels, and the
`scaled`/`weighted` flags. PIC requires a strictly bifurcating tree with finite
strictly positive branch lengths. Use `multi2di` explicitly for multifurcating
trees. PIC is not a general ancestral-state reconstruction API.

### Boundary Policies

- `pd` assigns zero PD to an empty community and the root-to-tip path length to
  a singleton when `include_root=true`.
- Empty or undefined community metrics return `NaN` rather than silently
  dropping rows; UniFrac and PhyloSor use the explicit diagonal rules above.
- Negative, infinite, or NaN community values raise `ArgumentError`.
- K/PIC reject missing and non-finite trait values.
- RaoD rejects non-ultrametric trees. PIC rejects multifurcations and zero or
  non-finite branch lengths.
- Randomized methods accept an injected `AbstractRNG`; seed-based `nri_nti`
  remains deterministic for the same seed.

## Performance Benchmarks

Command shape:

```julia
tree = read_tree(path)
distances_from_root(tree)
evol_distinct(tree)
sample = zeros(Int, 1, n)
sample[1, 1:div(n, 2)] .= 1
pd(sample, tree)
```

| Tips | read_tree | Root distances | ED | Single-site PD |
| ---: | ---: | ---: | ---: | ---: |
| 10,000 | 0.0170 s | 0.000045 s | 0.000209 s | 0.000807 s |
| 100,000 | 0.2786 s | 0.000302 s | 0.001825 s | 0.013961 s |
| 200,000 | 0.5099 s | 0.059488 s | 0.003887 s | 0.017693 s |
| 500,000 | 1.6090 s | 0.260128 s | 0.010796 s | 0.154415 s |

Measurements use balanced rooted trees with unit branch lengths on Julia 1.12.7
after a warm-up.

## Scope

EvoKits focuses on portable, tree-centered utilities:

- Newick I/O and `CompactTree` conversion
- Tree pruning, ultrametric conversion, and multifurcation resolution
- Tree identity and R/ape mapping
- ED, PD, SES-PD, SES-MPD, and SES-MNTD metrics
- Pagel's lambda and trait alignment
- Fritz-Purvis D for binary-trait phylogenetic signal
- EDGE2 scores and extinction-probability utilities
- Local-pool NRI/NTI with deterministic null sampling
