# EvoKits.jl

EvoKits is a Julia-native toolkit for practical phylogenetic work. It provides
tree I/O, a compact in-memory tree representation, topology utilities,
community phylogenetic metrics, tree identity mapping, and lightweight
comparative analysis.

It is an independent package extracted from the reusable tree-tooling layer of
EvoTraits. EvoKits does not include EvoTraits' continuous or discrete model
engines, ancestral reconstruction workflows, or their result types.

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

## Comparative Utility

```julia
trait = [1.2, 2.1, 1.8]
result = pagel_lambda(tree, trait)
summary_pagel_lambda(result)
```

`pagel_lambda` accepts a `CompactTree` and a univariate trait vector. Use
`align_traits_to_tree` when trait data must be reordered to match tree tips.

## Scope

EvoKits focuses on portable, tree-centered utilities:

- Newick I/O and `CompactTree` conversion
- Tree pruning, ultrametric conversion, and multifurcation resolution
- Tree identity and R/ape mapping
- ED, PD, SES-PD, SES-MPD, and SES-MNTD metrics
- Pagel's lambda and trait alignment

Continuous and discrete evolutionary-model fitting, reconstruction, SIMMAP,
and analysis workflows remain in EvoTraits.
