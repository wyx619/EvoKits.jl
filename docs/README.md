# EvoKits documentation

[Project home](../README.md) · [简体中文](README_zh.md)

Choose a guide by the analysis you want to perform. Each guide supplies its own
example inputs; run code blocks in order within that page.

## Start here

For community analysis, begin with [inputs and labels](community_input.md),
choose a [diversity or beta metric](community.md), then use
[regionalization and ordination](regionalization.md) to interpret a complete
distance matrix. The [project quick start](../README.md#quick-start) demonstrates
that sequence on a small dataset.

| Task | Guide |
| --- | --- |
| Read, query, prune or complete a tree; stream/sample posterior trees; slice time | [Trees](trees.md) |
| Import a matrix or long table; align species; export data; choose a reusable context | [Community inputs](community_input.md) |
| Calculate PD, endemism, diversity moments, Rao partition or pairwise/multiple-site beta | [Community diversity](community.md) |
| Compare average/nearest-taxon distances or branch overlap; request selected pairs | [Community distances](community_distances.md) |
| Choose a null model; calculate SES, beta SES or NRI/NTI | [Null models and community structure](null_models.md) |
| Summarize traits; aggregate species; form functional or lineage groups | [Traits and aggregation](aggregation.md) |
| Cluster sites; fit fixed-k GoM; average posterior distances; run/align NMDS | [Regionalization and ordination](regionalization.md) |
| Measure trait signal with lambda, K/PIC or D | [Phylogenetic signals](signals.md) |
| Calculate EDGE/EDGE2; select hotspots/coldspots; classify metric values | [Conservation metrics](conservation.md) |

NRI/NTI describe community structure relative to a species pool. Lambda, K and
D describe traits on a phylogeny; these analyses have separate guides.

## Reference and development

| Need | Reference |
| --- | --- |
| Ordering, output shapes, result fields and undefined values | [API conventions](api.md) |
| Choose a computation route; understand complexity, memory and benchmark conditions | [Performance](performance.md) |
| Find implementation files and their loading dependencies | [Source layout](development.md) |

Guides explain models and show usage. The API reference summarizes contracts;
the performance guide keeps algorithms, complexity and measurements together.
Use Julia's help mode, for example `?expected_pd`, for function docstrings.
