# Community diversity and structure

[Documentation](README.md) · [Home](../README.md) · [简体中文](community_zh.md)

On this page: [PD and bootstrap](#pd-variants-and-bootstrap) · [Diversity moments](#psv-psr-and-analytical-pd-moments) · [Rao](#rao-diversity-partition) · [Endemism and pairwise beta](#endemism-relative-pd-and-beta-decomposition) · [Multiple-site beta](#multiple-site-taxonomic-beta) · [Phylogenetic multiple-site beta](#multiple-site-phylogenetic-beta)

Choose a diversity or beta statistic, calculate it, and interpret its result.
The examples share this tree and labeled community:

```julia
using EvoKits, SparseArrays

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
community = [1 0 1 0; 0 1 1 0]
data = CommunityData(sparse(community);
    sample_labels = ["west", "east"], species_labels = tree.tip_labels)
```

## Community Metrics

| Question | Function | Section |
| --- | --- | --- |
| How much branch history does a site contain? | `pd`, `pd_calc`, `pd_bootstrap` | [PD variants](#pd-variants-and-bootstrap) |
| How variable, even or clustered are its species? | `psv`, `psr`, `pse`, `psc`, `psd` | [Diversity moments](#psv-psr-and-analytical-pd-moments) |
| What PD is expected under independent tip occurrence? | `expected_pd`, `pd_variance` | [Analytical PD moments](#psv-psr-and-analytical-pd-moments) |
| How is Rao diversity partitioned? | `raoD` | [Rao partition](#rao-diversity-partition) |
| Which sites contain restricted species or branches? | `weighted_endemism`, `phylo_endemism`, `relative_pd` | [Endemism and pairwise beta](#endemism-relative-pd-and-beta-decomposition) |
| How do sites differ? | `beta_diss`, `phylobeta` | [Pairwise beta](#endemism-relative-pd-and-beta-decomposition) |
| How diverse is a whole collection? | `beta_multi`, `phylobeta_multi` | [Multiple-site beta](#multiple-site-taxonomic-beta) |

### Export a long table and reconstruct the community

See [community inputs](community_input.md#export-a-long-table-and-reconstruct-the-community)
for constructors, sparse storage and lossless long-table round trips.

### Tree alignment and metric calls

See [alignment and reusable contexts](community_input.md#tree-alignment-and-metric-calls)
for missing species, label order and the choice of preparation function.

### PD variants and bootstrap

`pd_calc(community, tree; variant=:tbl)` returns `pd`, `richness`, and the
normalized variant, preserving site order. `pd_calc(tree, tips; variant=...)`
returns a scalar; `tips` are unique tip names or positions. Dense/CSC matrices
and `CommunityData` use positive abundances as presence.

| Variant | Original branches included |
| --- | --- |
| `:tbl` | Paths from observed tips to the original root. |
| `:mst` | Connections below the observed MRCA, excluding its common root stem. |
| `:ueh` | Branches with exactly one observed descendant. |
| `:sbl` | Branches with at least two observed descendants. |
| `:tip` | Original terminal branches of observed species. |

Empty sites give zero for every variant. A singleton has MST zero; UEH
includes its whole root path, while TIP includes only its original terminal edge.

```julia
using Random

variants = pd_calc(data, tree; variant = :ueh)
draws = pd_bootstrap(data, tree; reps = 99, rng = MersenneTwister(42))
draws.observed
draws.replicates              # sites × reps
fixed_richness = pd_bootstrap(tree, 2; reps = 99, method = :MST,
    tip_weights = [1.0, 2.0, 3.0, 4.0], rng = MersenneTwister(42))
```

Bootstrap samples without replacement from the complete tree-tip pool,
using optional weights and preserving each site's richness. Weights
follow tree-tip order, must be finite nonnegative Float64 values after conversion,
and must contain a positive value. Omitted weights are uniform. Tree-only
calls return a vector; community calls return `PDBootstrapResult`. Empty sites
give zero; `include_root=false` gives NaN for singleton community sites, as in `pd`.

If requested richness exceeds the number of positive-weight tips, remaining
slots are filled from zero-weight tips in original tip order. This differs
from R's rejection of insufficient support. Finite weights whose sum overflows
are scaled before normalization. R and Julia share the weighted sampling law
when positive support is sufficient, not identical random streams. See
[bootstrap costs](performance.md#pd-variants-and-bootstrap).

### PSV, PSR and analytical PD moments

`psv` measures the mean pair contrast `(Cii+Cjj)/2-Cij`. By default C is the
Brownian correlation matrix; `scale_vcv=false` uses covariance, so this contrast
is patristic distance / 2. `psr` multiplies PSV by site richness. Both functions
use tree summaries directly, including variance, and accept dense/CSC inputs
or `CommunityData`.

PSV variance describes uniform sampling **without replacement** at the
site's richness from all tree tips, including unobserved tips. It is zero at
full-pool richness. Unlike the leaf-order-dependent formula in the reference
picante implementation, it is invariant to tip order; raw-covariance variance
also includes unequal tip variances. Empty/singleton sites give NaN, and
`compute_var=false` leaves variance as NaN. Scaled calls require positive
root-to-tip lengths.

```julia
moment_data = CommunityData(sparse([1 1 0 0; 1 1 1 1]);
    species_labels=["A", "B", "C", "D"])
variability = psv(moment_data, tree)
variability.psv       # [1/2, 5/6]
variability.variance  # [1/18, 0]
richness_metric = psr(moment_data, tree)
moments = expected_pd(tree; upper_bound=false)
moments.expected_pd  # [1.875, 3.5, 4.875, 6.0]
moments.variance_pd  # [2.3671875, 2.375, 1.2421875, 0.0]
selected = expected_pd(tree; richness=[3, 1, 4], upper_bound=false)
selected.richness     # [3, 1, 4]
selected.expected_pd  # [4.875, 1.875, 6.0]
diversity = psd(moment_data, tree)
diversity.pse         # [1/2, 5/6]: equal abundances give PSE = PSV
diversity.psc         # [1/2, 1/2]
pse(sparse([2.0 1 0 0]), tree).pse  # [4/9]
```

PSV/PSR variance, including the corresponding PSD outputs, always uses all
original tree tips, including unobserved species. See [algorithm and storage
costs](performance.md#analytical-diversity-moments).

`pse` weights pair contrasts by relative abundance; `psc` averages each
present species' smallest `1-Cij` to another species. In raw covariance mode
PSC is one minus the mean largest shared branch length, can be negative,
and depends on branch-length units. In correlation mode, the best match on
a non-ultrametric tree can differ from the nearest tip by patristic distance.
Zero-length paths to other species remain eligible. `psd` returns all four
metrics, PSV/PSR variances and richness together, using the same definitions
and NaN rules as the individual calls. Dense/CSC and `CommunityData` inputs
are preserved. Algorithm and allocation details are in the [performance guide](performance.md#analytical-diversity-moments).

`expected_pd` and its alias `pd_variance` retain a different model: every tip
occurs **independently**, with probability `richness/tree.ntips`. Thus richness
is the expected number present. `upper_bound=false` gives the exact Bernoulli
variance; the default true gives picante's covariance upper bound. These are
not fixed-size rarefaction moments.

Omitting `richness` returns the full curve `1:tree.ntips`. Supply an integer
or a vector/range of integers in that interval to compute only the requested
points; their order and duplicates are retained. A scalar returns one-element
arrays; an empty vector returns empty arrays. Bool, noninteger types and
out-of-range values are rejected. Results own their arrays. Branch lengths
must be finite/nonnegative with a finite sum, including for an empty request.
See [algorithm and scaling](performance.md#analytical-diversity-moments).

### Rao diversity partition

`raoD` requires an ultrametric tree and uses patristic distance / 2, following
Picante. Its `RaoDResult` contains within-site `Dkk`, between-site `Dkl`,
`H[i,j] = Dkl[i,j] - (Dkk[i]+Dkk[j])/2`, pooled `total`, site-weighted `alpha`,
`beta = total-alpha`, and `Fst = beta/total`.

```julia
rao_data = CommunityData(sparse([2.0 1 0; 0 1 3]);
    sample_labels=["west", "east"], species_labels=["A", "B", "C"])
rao = raoD(rao_data, tree)
rao.Dkk           # [4/9, 3/4]
rao.total         # 8/7
rao.Fst           # 11/24
raoD(rao_data, tree; abundance_weighted=false)
rao_summary = raoD(rao_data, tree; pairwise=false)
rao_summary.Fst    # 11/24
size(rao_summary.Dkl), size(rao_summary.H) # ((0, 0), (0, 0))
```

The default `pairwise=true` returns both site-pair matrices. When only within-site
diversities and the pooled partition are needed, `pairwise=false` skips site-pair
calculations. It returns the same `RaoDResult` with independent `0×0` `Dkl` and
`H` matrices; `Dkk`, `total`, `alpha`, `beta` and `Fst` retain their definitions.

Abundance weights are the default. With `abundance_weighted=false`, species
within each site have equal weight, while pooled composition and site weights
still use original abundances. This existing Julia option can yield negative
beta/Fst; Picante has no corresponding option.

An empty site gives NaN Dkk and NaN Dkl/H rows and columns. It also makes
alpha/beta/Fst NaN, even with zero site weight; total can remain finite. All-empty
input gives NaN total. Zero input rows give empty arrays and alpha=0. Zero total
diversity gives undefined Fst (NaN when beta is zero). Dense/CSC matrices and
`CommunityData` share these rules and are left unchanged. See
[algorithm and storage costs](performance.md#rao-diversity-partition).

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

## Multiple-site taxonomic beta

`beta_multi` summarizes a selected collection with scalar `total`, `turnover`
and `nestedness`, using the multiple-site Sørensen or Jaccard definition.

```julia
collection = CommunityData(sparse([1 1 1 0; 1 0 0 1; 1 0 0 0]);
    sample_labels = ["s1", "s2", "s3"], species_labels = ["A", "B", "C", "D"])
overall = beta_multi(collection)
overall.total             # 3/5
overall.turnover          # 1/3
overall.nestedness        # 4/15
overall.richness          # [3, 2, 1]
overall.pooled_richness   # 4

subset = beta_multi(collection; samples = ["s3", "s1"], index_family = :jaccard)
subset.sample_labels      # ["s3", "s1"]
subset.sample_indices     # [3, 1], original input rows
subset.pooled_richness    # 3, counted within the selected collection
```

Dense/CSC matrices and `CommunityData` are accepted. All positive values denote
presence; abundance magnitudes do not weight the index. Matrix labels are set
through `sample_labels`, must be unique and nonempty, and default to
`sample_1`, ... . `samples` accepts a vector or tuple of distinct row indices
or string labels, preserving request order. Integer label identifiers become
strings; select those labels as strings. Bool masks and mixed index/label
selectors are rejected. All input values are validated, including excluded rows.
The result owns its labels, original row indices and selected-site richness.

With A equal to summed site richness minus pooled richness, and M/X the sums
of minimum/maximum unique counts across unordered site pairs, Sørensen total
is `(M+X)/(2A+M+X)` and turnover is `M/(A+M)`. Jaccard total is
`(M+X)/(A+M+X)` and turnover is `2M/(A+2M)`. Nestedness uses the equivalent
defining product to retain small values. These are collection-level indices;
averaging `beta_diss` entries gives a different statistic.

At least two selected sites are required. Empty sites remain in the collection;
all-empty or zero-column input yields three `NaN`s. With exactly one nonempty
site, total is one and turnover/nestedness are `NaN`. With no shared species
and at least two nonempty sites, total and turnover are one and nestedness zero.
Other cases follow their component denominators. No input is modified and no
full pairwise matrix is built; see the [algorithm notes](performance.md#multiple-site-taxonomic-beta).

## Multiple-site phylogenetic beta

`phylobeta_multi` applies the same multiple-site decomposition to branch
lengths. It accepts a dense/CSC matrix and tree, `CommunityData` and tree, or a
prepared `CommunityContext`. For the collection above:

```julia
phylogeny = mktempdir() do dir
    path = joinpath(dir, "phylogeny.nwk")
    write_tree(path, "((A:1,B:2):3,(C:4,D:5):6);")
    read_tree(path)
end
phylo_overall = phylobeta_multi(collection, phylogeny)
phylo_overall.total          # 17/31
phylo_overall.turnover       # 5/19
phylo_overall.nestedness     # 168/589
phylo_overall.pd             # [16, 15, 4]
phylo_overall.pooled_pd      # 21

phylo_context = prepare_community(collection, phylogeny)
phylo_subset = phylobeta_multi(phylo_context;
    samples = ["s3", "s1"], index_family = :jaccard)
phylo_subset.sample_indices  # [3, 1]
phylo_subset.pd              # [4, 16]
phylo_subset.pooled_pd       # 16, only branches occupied by selected sites
```

`samples` follows the selection and validation rules of `beta_multi`.
`PhylobetaMultiResult` owns the selected labels, original row indices, richness
and rooted PD; its three beta components are scalars. Positive abundances mean
presence. Matrix columns use tree-tip order or `tip_labels`; `CommunityData`
aligns species labels and zero-fills absent tree tips.

The formulas above use summed PD minus pooled PD for A, and unique branch
lengths for M/X. Pooled PD is the union of branches occupied by the selected
sites, including paths to the original root. Unobserved species contribute no
branches. This agrees with betapart's direct tree calculation. The local
phyloregion core uses full tree length for pooled PD, which can give negative
multi-site components when branches are unobserved; that behavior is not used.

Context calls reuse branch sets and recompute occupancy for the selected sites.
They do not change the full-collection occupancy used by WE/PE. At least two
selected sites are required. Empty or zero-PD sites remain; all-zero PD gives
three `NaN`s. Exactly one positive-PD site gives total one and `NaN`
turnover/nestedness. Other cases follow the component denominators. With two
sites the decomposition reduces to `phylobeta`. Inputs remain unchanged and no
full pairwise matrix is built; see the [algorithm notes](performance.md#multiple-site-phylogenetic-beta).


### Distances between communities

See [community distances](community_distances.md) for `comdist`, `comdistnt`,
UniFrac/PhyloSor, selected pairs and repeated distance calculations.

## NRI and NTI

See [null models and community structure](null_models.md#nri-and-nti) for
local species pools, sign conventions and batched NRI/NTI.

## Phylogenetic beta SES

See [beta SES](null_models.md#phylogenetic-beta-ses) for randomization models,
valid-draw counts and the interpretation of z scores and p-values.
