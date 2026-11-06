# Source layout

[Documentation](README.md) · [Home](../README.md) · [简体中文](development_zh.md)

EvoKits has one Julia module. [EvoKits.jl](../src/EvoKits.jl) groups public
exports by function and loads source files in dependency order. Subdirectories
organize responsibilities within that module; they do not introduce submodules.

## Find an implementation

| Directory under `src/` | Responsibility | Entry point |
| --- | --- | --- |
| `base/` | Shared labels, tree representation, I/O and numerical helpers | [labels.jl](../src/base/labels.jl), [tree_io.jl](../src/base/tree_io.jl) |
| `tree/queries/` | Validation, queries, distances and clade groups | [clade_groups.jl](../src/tree/queries/clade_groups.jl) |
| `tree/operations/` | Pruning, transformations, time slices, completion and sampling | [phylobuilder.jl](../src/tree/operations/phylobuilder.jl) |
| `tree/mapping/` | Node/edge identities and R/ape numbering | [mapping.jl](../src/tree/mapping/mapping.jl) |
| `tree/metrics/` | Evolutionary distinctiveness and tree shape | [metrics.jl](../src/tree/metrics/metrics.jl) |
| `traits/` | Trait alignment and phylogenetic signal | [alignment.jl](../src/traits/alignment.jl) |
| `community/core/` | Inputs, long-table export, branch summaries, contexts and adapters | [input.jl](../src/community/core/input.jl), [data_metrics.jl](../src/community/core/data_metrics.jl) |
| `community/diversity/` | PD, species diversity, distances, endemism, Rao and beta | [community_distances.jl](../src/community/diversity/community_distances.jl), [phylogenetic_species.jl](../src/community/diversity/phylogenetic_species.jl) |
| `community/aggregation/` | Trait/species/lineage groups and functional beta | [functional_beta.jl](../src/community/aggregation/functional_beta.jl) |
| `community/regionalization/` | Clustering, GoM, posterior means, NMDS and indicators | [regionalization.jl](../src/community/regionalization/regionalization.jl) |
| `community/null_models/` | Randomization, SES, beta SES and NRI/NTI | [ses.jl](../src/community/null_models/ses.jl) |
| `community/conservation/` | EDGE/EDGE2, hotspots/coldspots and metric classes | [conservation.jl](../src/community/conservation/conservation.jl) |

## Loading and shared machinery

MKL loads before other dependencies. Base and tree definitions precede traits
and communities. Shared label validation lives in `base/labels.jl`, allowing
tree completion to work independently of community input. Load phases can span
directories; files in a directory need not be loaded consecutively.

| Work shared across features | Files |
| --- | --- |
| Stream framing, NEXUS commands and selected-tree sampling | [tree_stream.jl](../src/base/tree_stream.jl), [tree_nexus.jl](../src/base/tree_nexus.jl), [tree_sampling.jl](../src/tree/operations/tree_sampling.jl) |
| Presence-based branch sets and occupancy for PD/PE/beta | [branches.jl](../src/community/core/branches.jl), [context.jl](../src/community/core/context.jl) |
| Sparse records, automatic connecting trees and distance reuse | [community_distances.jl](../src/community/diversity/community_distances.jl), [comdist_virtual.jl](../src/community/diversity/comdist_virtual.jl), [distance_context.jl](../src/community/core/distance_context.jl) |
| Site trees for PSV/PSE/PSC and PD variants | [species_virtual.jl](../src/community/diversity/species_virtual.jl), [pd_workspace.jl](../src/community/diversity/pd_workspace.jl) |
| Existing-distance and direct posterior-beta accumulation | [mean_dist.jl](../src/community/regionalization/mean_dist.jl), [mean_beta.jl](../src/community/regionalization/mean_beta.jl) |

Community adapters load after the kernels they call. Rao reuses the distance
records/topology; PSD shares site preparation; direct posterior beta reuses mean
accumulators but keeps per-tree buffers separate from single-tree contexts.
See [performance](performance.md) for algorithms and complexity, and [API
conventions](api.md) for snapshot/ordering contracts.

## Make and validate a change

Keep numerical tests in `test/` and small offline R reference tables in
`test/fixtures/`. Fixture READMEs describe their provenance; ordinary tests
should not require an R installation or network access.

Run affected module tests for behavior changes. For source reorganization,
also check package loading and dependent modules. For documentation edits,
check relative links, bilingual parity and changed runnable examples; a full
numerical suite is not required for documentation-only changes. Update the
task guide for usage, `api.md` for contracts and `performance.md` for costs.

Local validation scripts, historical snapshots and detailed acceptance logs
live in the ignored `validation/` directory. They are not published API or
documentation dependencies. Historical Git paths refer to their original commit.
