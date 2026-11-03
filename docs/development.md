# Source layout

[Home](../README.md) · [简体中文](development_zh.md)

EvoKits uses one Julia module. [EvoKits.jl](../src/EvoKits.jl) declares the
public exports and loads source files in dependency order. Subdirectories group
related code within that module; exported names and calling conventions stay
the same.

| Directory | Responsibility | Representative source |
|---|---|---|
| `base/` | Shared identifiers, tree representation, I/O and numeric helpers | [labels.jl](../src/base/labels.jl) |
| `tree/queries/` | Validation, queries, distances and clade grouping | [clade_groups.jl](../src/tree/queries/clade_groups.jl) |
| `tree/operations/` | Pruning, transformation, multifurcations, time slices and completion | [phylobuilder.jl](../src/tree/operations/phylobuilder.jl) |
| `tree/mapping/` | Node/edge identities and R/ape numbering | [mapping.jl](../src/tree/mapping/mapping.jl) |
| `tree/metrics/` | Evolutionary distinctiveness and tree shape | [metrics.jl](../src/tree/metrics/metrics.jl) |
| `traits/` | Trait alignment and phylogenetic signal | [alignment.jl](../src/traits/alignment.jl) |
| `community/core/` | Community input, branch traversal, contexts and metric adapters | [input.jl](../src/community/core/input.jl) |
| `community/diversity/` | PD, MPD/MNTD, endemism and beta diversity | [beta.jl](../src/community/diversity/beta.jl) |
| `community/aggregation/` | Species/lineage aggregation, functional grouping and count selection | [functional_beta.jl](../src/community/aggregation/functional_beta.jl) |
| `community/regionalization/` | Partitions, Infomap, linkage evaluation, distance means, NMDS and indicators | [regionalization.jl](../src/community/regionalization/regionalization.jl) |
| `community/null_models/` | Randomization, SES, beta SES and NRI/NTI | [ses.jl](../src/community/null_models/ses.jl) |
| `community/conservation/` | EDGE/EDGE2, hotspots/coldspots and value classification | [conservation.jl](../src/community/conservation/conservation.jl) |

Shared label validation lives in `base/labels.jl`: nonempty strings or integer
identifiers become unique string labels, with the existing default-label and
error rules. Tree completion and community input call the same helpers. The
base and tree layers load before traits and community; `phylobuilder` therefore
works without community definitions. Dependencies between community files
determine their load order even when they occupy different subdirectories.

`infomap.jl` handles community validation and site results; `infomap_backend.jl`
isolates IGraphs native objects, random state and calls. EvoKits directly depends
on IGraphs; Graphs is its indirect dependency.

Keep numerical tests under `test/` and small offline reference tables under
`test/fixtures/`. Run the affected test files and package-loading checks for
structural changes. Local validation scripts and historical source snapshots
live in the ignored `validation/` directory; Git paths inside historical
`git show` commands refer to the layout of their original commit.
