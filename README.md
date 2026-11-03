# EvoKits.jl

[简体中文](README_zh.md)

Julia tools for phylogenetic trees, community diversity, biogeographic
regionalization and comparative analysis, with dense and sparse community inputs.

## Install

```julia
using Pkg
Pkg.add(url = "https://github.com/wyx619/EvoKits.jl")
```

## Quick start

Create a small tree and a labeled community, then calculate phylogenetic diversity
(PD), phylogenetic endemism (PE) and a two-region partition:

```julia
using EvoKits, SparseArrays

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,C:2);")
    read_tree(path)
end

data = CommunityData(
    sparse([1 0 1; 0 1 1; 1 1 0; 1 1 1]);
    sample_labels = ["north", "south", "east", "west"],
    species_labels = ["A", "B", "C"],
)

pd(data, tree).pd
phylo_endemism(data, tree).pe

beta = phylobeta(data, tree)
regions = phyloregion(beta; k = 2)
regions.membership
```

Rows are sites and columns are species. `CommunityData` matches species to tree
tips by label.

## Guides

| Topic | Contents |
| --- | --- |
| [Trees](docs/trees.md) | Newick I/O, topology, pruning, R/ape mapping and branch time slices |
| [Community diversity](docs/community.md) | Sparse inputs, PD, WE/PE, relative PD, beta decomposition, SES and reusable contexts |
| [Traits and aggregation](docs/aggregation.md) | Trait summaries, species groups, k-prototypes, functional beta and lineage aggregation |
| [Regionalization and ordination](docs/regionalization.md) | Eight linkage methods, partition evaluation, posterior-tree mean distances and NMDS |
| [Phylogenetic signals](docs/signals.md) | NRI/NTI, Pagel's lambda, K, independent contrasts and D |
| [Conservation metrics](docs/conservation.md) | EDGE2 scores and extinction probability tools |
| [API conventions](docs/api.md) | Combined example, ordering, result types and boundary behavior |
| [Performance](docs/performance.md) | Runtime backend and benchmark timings |
