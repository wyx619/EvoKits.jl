# Source layout

[Home](../README.md) · [简体中文](development_zh.md)

EvoKits uses one Julia module. [EvoKits.jl](../src/EvoKits.jl) declares the
public exports and loads source files in dependency order. Subdirectories group
related code within that module; exported names and calling conventions stay
the same.

The entry file groups exports by function, with result types next to their
functions. NRI/NTI and SES belong to null models and community structure;
phylogenetic signal and conservation scores have separate groups. Source-loading
sections describe dependency phases, which may span several directories. Keep
MKL as the first loaded dependency.

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
| `community/regionalization/` | Partitions, linkage evaluation, distance means, NMDS and indicators | [regionalization.jl](../src/community/regionalization/regionalization.jl) |
| `community/null_models/` | Randomization, SES, beta SES and NRI/NTI | [ses.jl](../src/community/null_models/ses.jl) |
| `community/conservation/` | EDGE/EDGE2, hotspots/coldspots and value classification | [conservation.jl](../src/community/conservation/conservation.jl) |

Shared label validation lives in `base/labels.jl`: nonempty strings or integer
identifiers become unique string labels, with the existing default-label and
error rules. Tree completion and community input call the same helpers. The
base and tree layers load before traits and community; `phylobuilder` therefore
works without community definitions. Dependencies between community files
determine their load order even when they occupy different subdirectories.

`base/tree_stream.jl` frames and parses multi-tree Newick records;
`tree/operations/tree_sampling.jl` samples them in two scans, parsing only
selected records. Both load before traits and community and have no dependency
on the community layer.

`gom.jl` exposes fixed-k mixture results and the indicator adapter;
`gom_fit.jl` implements sparse likelihood/EM accumulation and safeguarded
log-ratio acceleration. It uses existing Julia dependencies.

`nmds_alignment.jl` matches labels and fits orthogonal Procrustes transforms;
`nmds_colors.jl` maps coordinates to the fixed phyloregion palette and keeps
reusable bounds. Both use existing dependencies and operate on copied scores.

Keep numerical tests under `test/` and small offline reference tables under
`test/fixtures/`. Run the affected test files and package-loading checks for
structural changes. Local validation scripts and historical source snapshots
live in the ignored `validation/` directory; Git paths inside historical
`git show` commands refer to the layout of their original commit.
