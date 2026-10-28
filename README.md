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

For Caper-compatible ED handling on polytomies, use
`evol_distinct(tree; type = :fair_proportion,
polytomy_correction = :isaac)`. The other supported modes are `:mooers` and
`:none`.

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
- `ses_pd`, `ses_mpd`, and `ses_mntd` return typed results containing observed
  values, null means/SDs, z scores, rank p-values, run count, and null model.
  Empty or singleton MPD/MNTD observations are `NaN`.
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
