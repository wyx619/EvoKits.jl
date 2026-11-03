# Community diversity and endemism

[Home](../README.md) · [简体中文](community_zh.md)

Labeled dense and sparse inputs, PD, endemism, beta decomposition, reusable contexts and null models.

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
