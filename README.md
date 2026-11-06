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

Browse the [documentation index](docs/README.md) for inputs, distances, null models and task-based reading paths.

| Topic | Contents |
| --- | --- |
| [Trees](docs/trees.md) | Newick I/O, topology, pruning, taxonomy-aided completion, R/ape mapping and time slices |
| [Community diversity and structure](docs/community.md) | PD, WE/PE, relative PD, diversity moments, Rao partition and beta decomposition |
| [Traits and aggregation](docs/aggregation.md) | Trait summaries, k-prototypes, group-count selection, functional beta and lineage aggregation |
| [Regionalization and ordination](docs/regionalization.md) | Eight linkage methods, mixed memberships, posterior means, NMDS alignment/colors and KL indicators |
| [Phylogenetic signals](docs/signals.md) | Pagel's lambda, K, independent contrasts and D |
| [Conservation metrics](docs/conservation.md) | EDGE/EDGE2, extinction probabilities, hotspots/coldspots and Jenks classes |
| [API conventions](docs/api.md) | Ordering, result shapes/fields and boundary behavior |
| [Performance](docs/performance.md) | Computation routes, complexity, memory and benchmark conditions |
| [Source layout](docs/development.md) | Internal directories and dependency order |
