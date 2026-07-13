# EvoKits.jl

EvoKits is a Julia-native toolkit for practical phylogenetic work. It provides
tree I/O, a compact in-memory tree representation, topology utilities,
community phylogenetic metrics, tree identity mapping, and lightweight
comparative analysis.

It is an independent package extracted from the reusable tree-tooling layer of
EvoTraits. EvoKits does not include EvoTraits' continuous or discrete model
engines, ancestral reconstruction workflows, or their result types.

The stable public conventions are documented in [API.md](API.md): tip and
community ordering, result fields, RNG behavior, and empty/singleton policies.

## Install and Load

For local development:

```julia
julia --project=.

julia> using Pkg
julia> Pkg.instantiate()
julia> using EvoKits
```

EvoKits loads MKL as its first dependency. This selects Intel MKL as Julia's
BLAS backend for the current Julia session.

## Tree I/O and Topology

```julia
using EvoKits

raw_tree = load_newick_tree("example.nwk")
tree = to_compact_tree(raw_tree)

ultrametric_tree = make_ultrametric(tree)
binary_tree = multi2di(ultrametric_tree; random = false)
save_newick_tree("binary.nwk", binary_tree)
```

`CompactTree` is the common input type for EvoKits tools. It caches topology,
edge lengths, traversal order, and tip-label lookup while preserving a clear
Newick conversion boundary through `to_compact_tree` and `from_compact_tree`.

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

## Comparative Utility

```julia
trait = [1.2, 2.1, 1.8]
result = pagel_lambda(tree, trait)
summary_pagel_lambda(result)
```

`pagel_lambda` accepts a `CompactTree` and a univariate trait vector. Use
`align_traits_to_tree` when trait data must be reordered to match tree tips.

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

save_newick_tree("example.nwk", "((A:1,B:1):1,(C:1,D:1):1);")
tree = to_compact_tree(load_newick_tree("example.nwk"))

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

Continuous and discrete evolutionary-model fitting, reconstruction, SIMMAP,
and analysis workflows remain in EvoTraits.

Reference checks and scale smoke benchmarks are recorded in
`PERFORMANCE_BENCHMARKS.md`; the current test suite includes Picante/Caper
hand fixtures and explicit empty/singleton/invalid-tree policies.
