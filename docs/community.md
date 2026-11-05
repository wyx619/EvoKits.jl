# Community diversity and structure

[Home](../README.md) · [简体中文](community_zh.md)

Labeled dense and sparse inputs, PD, endemism, beta decomposition, NRI/NTI, reusable contexts and null models.

The examples on this page share the following tree:

```julia
using EvoKits

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
```

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

### Export a long table and reconstruct the community

`community_longtable(data)` returns a named tuple with a `table` (`DataFrame`)
and complete `sample_labels` / `species_labels` vectors. Each positive
site–species combination has one row, with fixed `sample`, `species`, and
`value` columns. Records follow species-column order, then site-row order.
Abundances and fractional weights retain their values and numeric type;
integer counts stay in the value column. Zero entries are omitted.

```julia
using SparseArrays

export_source = CommunityData(sparse([0.0 2.5 0; 0 0 0; 3 0.25 0]);
    sample_labels=["s2", "empty", "s1"], species_labels=["C", "A", "absent"])
records = community_longtable(export_source)
records.table
restored = CommunityData(records.table;
    sample_labels=records.sample_labels, species_labels=records.species_labels)
restored.matrix == export_source.matrix # true
restored.sample_labels                 # ["s2", "empty", "s1"]
restored.species_labels                # ["C", "A", "absent"]
```

Keep both label vectors with the table when saving it: nonzero records alone
lose empty sites, all-zero species columns, and possibly site order. This also
applies to an entirely empty table or a matrix with no rows/columns. Exported
arrays are independent of the input. Reconstruction uses CSC storage and the
existing long-table rules, including conversion of Bool values to integer
counts. CSC export reads stored entries directly; see [storage costs](performance.md#community-long-table-export).

### Tree alignment and metric calls

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
    1 0 1 0
    0 1 1 0
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

Rows are communities and columns are tips. Pass `tip_labels` when input
columns are not already in tree-tip order. Empty communities return `NaN`
where a distance or diversity value is undefined; singleton communities are
accepted by PD, UniFrac, PhyloSor, and RaoD.

### PSV, PSR and analytical PD moments

`psv` measures the mean pair contrast `(Cii+Cjj)/2-Cij`. By default C is the
Brownian correlation matrix; `scale_vcv=false` uses covariance, so this contrast
is patristic distance / 2. `psr` multiplies PSV by site richness. Both functions
use tree summaries directly, including variance, and accept dense/CSC inputs
or `CommunityData`.

PSV variance now describes uniform sampling **without replacement** at the
site's richness from all tree tips, including unobserved tips. It is zero at
full-pool richness. This corrects the previous leaf-order-dependent formula
shared with picante; raw-covariance variance also includes unequal tip
variances. Empty/singleton sites give NaN, and `compute_var=false` leaves
variance as NaN. Scaled calls require positive root-to-tip lengths.

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
diversity = psd(moment_data, tree)
diversity.pse         # [1/2, 5/6]: equal abundances give PSE = PSV
diversity.psc         # [1/2, 1/2]
pse(sparse([2.0 1 0 0]), tree).pse  # [4/9]
```

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
not fixed-size rarefaction moments. Branch lengths must be finite/nonnegative
with a finite sum. See [algorithm and scaling](performance.md#analytical-diversity-moments).

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
```

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

### Distances between communities

`comdist` averages all cross-community tip distances. Its diagonal includes
zero self-distances as well as distances between different tips. Positive
values supply equal species weights by default; `abundance_weighted=true`
uses their abundances. Empty sites give `NaN` rows and columns.

`comdistnt` averages the nearest-target distance in each direction. Directions
receive equal weight by default; with `abundance_weighted=true`, source
abundances weight each directional mean and original site totals weight the
two means. A source tip without an eligible target is omitted and remaining
weights are renormalized. Use the valid direction if only one exists, or
`NaN` if neither does. Nonempty diagonals are zero by default.

```julia
distance_data = CommunityData(sparse([2.0 1 0; 0 1 3; 0 0 0]);
    sample_labels=["west", "east", "empty"], species_labels=["A", "B", "C"])
cross_mean = comdist(distance_data, tree; abundance_weighted=true)
nearest = comdistnt(distance_data, tree; exclude_conspecifics=true)
cross_mean[3, :] # [NaN, NaN, NaN]
nearest[1, 1]   # 2.0: nearest different species within the first site
requested_pairs = [(1, 2), (2, 1), (1, 1), (3, 2)]
selected_mean = comdist(distance_data, tree;
    abundance_weighted=true, pairs=requested_pairs)
# [10/3, 10/3, 8/9, NaN]
selected_nearest = comdistnt(distance_data, tree;
    abundance_weighted=true, exclude_conspecifics=true, pairs=requested_pairs)
# [22/7, 22/7, 2, NaN]
```

Both `comdist` and `comdistnt` accept `pairs` and return a vector in request order, including
repeated, reversed and self pairs. Pair indices refer to the complete input
site order; an empty request returns an empty vector. Values, labels and tree
lengths are validated even for sites absent from the requests. Only requested
endpoints are retained for computation, and repeated/reversed requests reuse
one evaluation. `comdistnt` computes both directional means for each unique
pair, reusing each needed target pass. To average off-diagonal pairs with `mean_dist`, supply labeled
tuples containing `distances`, `sample_labels` and `pairs`, and explicit
output `pairs`. `mean_dist` requires self-distances to be zero or NaN, while
`comdist` self-pairs include within-site diversity. `comdistnt` self-pairs are
zero by default and may be nonzero with `exclude_conspecifics=true`.

`exclude_conspecifics=true` excludes species identity. Different species can
still have distance zero on zero-length branches; singleton diagonals then
give `NaN`. These are the existing Julia conventions. Picante combines
directional means differently when richness or site totals differ and excludes
all zero distances; it also returns a `dist` object with zero diagonals.

Both functions accept dense/CSC matrices or `CommunityData`, preserve site
order, and require finite non-negative branch lengths with a finite sum.
They leave inputs unchanged and support column alignment through `tip_labels`
or `CommunityData` labels. Their CSC paths avoid a dense community and a full
tip-distance matrix. Default output is a full site-pair matrix; selected
distance output contains only the requested values. See
[algorithm and storage costs](performance.md#community-phylogenetic-distances).

For repeated distance calculations, prepare a shared snapshot:

```julia
distance_context = prepare_community_distances(distance_data, tree)
reused_mean = comdist(distance_context; abundance_weighted=true,
    pairs=requested_pairs)
reused_nearest = comdistnt(distance_context; abundance_weighted=true,
    exclude_conspecifics=true, pairs=requested_pairs)

fixed_context = prepare_community_distances(distance_data, tree;
    pairs=requested_pairs)
fixed_nearest = comdistnt(fixed_context; abundance_weighted=true,
    exclude_conspecifics=true) # same request order as selected_nearest
```

`CommunityDistanceContext` stores both abundance and presence weights, one
owned tree topology and the workspaces for both distances. With no preparation
`pairs`, later calls can request full matrices or change the pairs; these calls
reuse the topology for all prepared sites. Fixing `pairs` during preparation
retains only their nonempty endpoints and caches request schedules. Such calls
default to the fixed vector request; rebuild to change its scope or order.

All input values, species labels, branch lengths and site labels are checked
at preparation; site labels must be unique and nonempty. The original matrix
and tree are not retained, and later edits do not affect the snapshot. Rebuild
for changed data. Treat context fields as internal and read-only. Results own
their arrays; calls sharing a context are serialized to protect its reusable
workspaces. No distance matrix is cached. The presence-based `CommunityContext`
for PD/PE/beta retains different summaries and remains a separate entry point.

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

## NRI and NTI

NRI and NTI describe within-community phylogenetic structure relative to a
null model. They are the negative standardized effects of MPD and MNTD:
`NRI = (null_mpd_mean - observed_mpd) / null_mpd_sd`, with the corresponding
MNTD formula for NTI. Positive values indicate clustering and negative values
indicate overdispersion relative to the chosen null model.

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
derived per grid, independent of processing order:

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
in parallel when Julia starts with multiple threads. Replace the example
paths below with your own input files.

```julia
batch = nri_nti_table(tree, "local_pool.tsv", "observed_community.tsv";
                        nreps = 999, seed = 20260713, threaded = true,
                        delim = '\t')
batch.grid_ids
batch.results
```

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

## Phylogenetic beta SES

`phylobeta_ses` standardizes one phylogenetic beta component against a null
model. It reuses the sparse randomization and online statistics of the site
SES methods. Choose `:sorensen` or `:jaccard` and one of `:total`, `:turnover`,
or `:nestedness` per call.

```julia
using Random

effect = phylobeta_ses(comm, tree; index_family = :sorensen,
    component = :turnover, null_model = :taxa_labels, runs = 999,
    pairs = [(1, 2)], rng = MersenneTwister(42))
effect.observed
effect.null_mean
effect.null_sd
effect.z_score
effect.p_value
effect.valid_replicates
```

Dense/CSC matrices and `CommunityData` are accepted; positive abundances
denote presence. The complete tree-tip pool remains eligible, including
unobserved species. `CommunityData` aligns by species name and fills missing
tree species with zeros; matrix inputs require a column for every tree tip.

| Null model | Randomization and constraints |
| --- | --- |
| `:taxa_labels` (`:tipshuffle`) | One shared permutation of tree-tip labels; preserves all taxonomic pairwise overlaps. |
| `:richness`, `:sample_pool` | Independently samples each site's richness from the complete tree-tip pool. |
| `:phylogeny_pool` | Richness sampling followed by a shared label permutation. |
| `:frequency` | Preserves each species' site occupancy while randomizing its sites; site richness can change. |
| `:trialswap`, `:independentswap` | Checkerboard trials preserve both margins and restart from the observed community for each replicate. |

Swap effort is controlled by `iterations` and `max_attempts_per_swap`; finite
swap runs do not guarantee uniform sampling. Use fresh RNGs with the same
seed and settings to compare components under identical draws.

`PhylobetaSESResult` contains the six arrays shown above, site labels and
richness, requested pairs, and null-model settings. Arrays are symmetric full
matrices by default. `pairs` returns vectors in request order, retaining
repeated, reversed, and self pairs. Every draw randomizes the complete input;
only branch sets for requested endpoints are built. Replicate matrices are
never retained, and selected-pair output avoids full site-by-site matrices.

Means and sample SDs use only finite draws, counted per pair. No valid draw
gives `NaN` mean/SD/z/p; one valid draw gives `NaN` SD/z. Constant nulls retain
SD zero and yield `NaN` z. Nonfinite observed beta gives `NaN` z/p even when
finite null summaries exist. Empty and zero-PD sites follow `phylobeta`,
including undefined diagonals.

The lower-tail mid-rank score is
`p = (1 + count(null < observed) + 0.5count(null == observed)) / (valid + 1)`.
Ties use exact computed values. This is not a two-sided p-value; its
denominator uses valid draws, while older site SES retains `runs + 1`.
R's `rowwise` and `colwise` label reorderings are undone by alignment, so
these models are rejected. R's replacement of zero SD by one is not used.
SES z scores can be negative or undefined and are unsuitable as distances
for `phyloregion`.
