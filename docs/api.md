# API conventions

[Documentation](README.md) · [Home](../README.md) · [简体中文](api_zh.md)

On this page: [Ordering](#shared-ordering) · [Matrices](#matrix-semantics) · [Community results](#community-result-types) · [Other results](#grouping-and-conservation-results) · [Trait results](#trait-result-types) · [Boundaries](#boundary-policies)

Use this page to check ordering, result shapes and boundary behavior. Models
and examples are in the [task guides](README.md); algorithm costs and measured
performance are in the [performance guide](performance.md).

## End-to-End Example

The [quick start](../README.md#quick-start) combines labeled community input,
PD, PE and regionalization. See [community inputs](community_input.md) for
constructors and [phylogenetic signals](signals.md) for trait analysis.

## Public API Contract

### Shared Ordering

| Input or result | Ordering |
| --- | --- |
| Tree | `tree.tip_labels` is the canonical tip order. Tip positions are `1:tree.ntips`; compact node IDs are `1:tree.nnodes`. |
| Community matrix | Rows are sites, columns are species. Without `tip_labels`, columns follow tree-tip order. With it, labels describe input columns and those columns are aligned to the tree. Direct matrix calls require all tree tips. |
| `CommunityData` | Species align by name; omitted tree species receive zero columns. Unknown species fail, including all-zero columns. Site order is preserved. |
| Trait vector | Follows tree-tip order unless `tip_labels` describes a different input order. |
| Site result | Vectors and matrix axes retain input site order unless an explicit subset is requested. |

Community values must be finite nonnegative real numbers. Positive values
denote presence; abundance-weighted functions use their magnitudes. Labels
must be unique and nonempty. Integer identifiers become strings. Tip-label
lookup uses `tree.tipname_to_id`; a node ID is not a tip position.

Constructors and preparation functions own their arrays. Treat context fields
as read-only snapshots; rebuild after changing inputs. Returned metric arrays
can be modified independently. See [context selection](community_input.md#repeated-analyses)
and the [distance-context contract](community_distances.md#distances-between-communities).

### Matrix Semantics

| Function | Default output and diagonal | Selected `pairs` |
| --- | --- | --- |
| `comdist` | `n×n` mean cross-species distances; diagonal includes zero self-distances and is generally nonzero | Vector in request order |
| `comdistnt` | `n×n` symmetric nearest-taxon means; nonempty diagonal zero by default; exclusion gives within-site nearest-other-species means | Vector in request order |
| `unifrac` | `n×n` branch dissimilarities; diagonal zero | Use `phylobeta` for selected beta pairs |
| `phylosor` | `n×n` branch similarities; diagonal one for positive PD, zero for zero PD | Use `phylobeta` for selected beta pairs |
| `beta_diss`, `phylobeta` | Three symmetric `n×n` component arrays | Three vectors in request order |
| `phylobeta_ses` | Symmetric observed/null/SES arrays | Vectors in request order |
| `mean_dist` | Symmetric distances and finite-observation counts | Vectors in request order |

`pairs` use complete input site indices; repeated, reversed and self pairs
are retained. Empty requests return empty vectors and still validate inputs.
Matrix-returning calls use ordinary `Float64` arrays; count arrays are integers.
See [distance definitions](community_distances.md), [beta decomposition](community.md#endemism-relative-pd-and-beta-decomposition)
and [null models](null_models.md).

`community_longtable(data)` returns `(table, sample_labels, species_labels)`.
The table contains `sample`, `species`, `value`, with one row per positive
cell in species-then-site order. Preserve both complete label lists to retain
empty sites and all-zero species on [reconstruction](community_input.md#export-a-long-table-and-reconstruct-the-community).

### Community Result Types

| Function | Result | Main fields |
| --- | --- | --- |
| `pd` | `PDResult` | `sample_labels`, `pd`, `richness`, `include_root` |
| `weighted_endemism` | `WeightedEndemismResult` | `sample_labels`, `we`, `richness`, `species_occupancy` (input species order) |
| `phylo_endemism` | `PhyloEndemismResult` | `sample_labels`, `pe`, `richness`, `branch_occupancy` (tree edge order), `weighted` |
| `relative_pd` | `RelativePDResult` | `sample_labels`, `rpd`, `pd_normalized`, `pd_equal_normalized`, `richness`, `total_branch_length`, `positive_branch_count` |
| `psv` | `PSVResult` | `psv`, `richness`, `variance`, `scale_vcv` |
| `psr`, `pse`, `psc` | Named tuples | Metric array and `richness`; PSR also returns `variance` |
| `psd` | Named tuple | `psv`, `var_psv`, `psc`, `psr`, `var_psr`, `pse`, `richness` |
| `expected_pd`, `pd_variance` | `ExpectedPDResult` | `richness`, `expected_pd`, `variance_pd`, `upper_bound` |
| `raoD` | `RaoDResult` | Within-site `Dkk`, pairwise `Dkl`/`H`, pooled `total`, `alpha`, `beta`, `Fst` |
| `pd_calc` | Scalar for tree/tip subset; named tuple for communities | `pd`, `richness`, normalized `variant` |
| `pd_bootstrap` | Vector for tree/richness; `PDBootstrapResult` for communities | `observed`, sites × reps `replicates`, `richness` |
| `beta_diss`, `phylobeta` | `BetaDissResult`, `PhylobetaResult` | `sample_labels`, `total`, `turnover`, `nestedness`, `index_family`, `pairs` |
| `beta_multi` | `BetaMultiResult` | Scalar components, `index_family`, selected `sample_labels`, original `sample_indices`, `richness`, `pooled_richness` |
| `phylobeta_multi` | `PhylobetaMultiResult` | Scalar components, `index_family`, selected `sample_labels`, original `sample_indices`, `richness`, rooted `pd`, selected-union `pooled_pd` |
| `ses_pd`, `ses_mpd`, `ses_mntd` | Typed site SES results | `sample_labels`, `ntaxa`, prefixed observed/null/rank/z/p arrays and run settings; MPD/MNTD also supply `mpd_nri`/`mntd_nti` |
| `phylobeta_ses` | `PhylobetaSESResult` | `sample_labels`, `richness`, `observed`, `null_mean`, `null_sd`, `z_score`, `p_value`, `valid_replicates`, `pairs` and run settings |

`expected_pd` defaults to `1:tree.ntips`; `richness` accepts an integer or an
integer vector/range within that interval. All three result arrays preserve
request order and duplicates, including one-element and empty requests.
Its independent-occurrence model differs from the fixed-size sampling model
of PSV variance; see [diversity moments](community.md#psv-psr-and-analytical-pd-moments).

`raoD(...; pairwise=false)` returns independent `0×0` `Dkl`/`H` matrices;
`Dkk` and the four scalars keep their definitions. See [Rao partition rules](community.md#rao-diversity-partition).
Multiple-site beta returns an overall scalar decomposition, not an average
of pairwise beta.

### Posterior tree input

`read_trees(paths; burnin=0, thin=1)` returns a single-pass `TreeStream` of
`CompactTree`s. Use its do-block to close files even on early exit. Filtering
and Newick/NEXUS detection are per file. NEXUS record numbers count only
`TREE`/`UTREE` definitions; translations apply within their `TREES` block.

`tree_sampler` returns iterable `TreeSampleResult`: independent `trees` and
draw-order `sample_indices`, `source_paths`, `source_indices`, `source_lines`,
plus `population_size`, `replace`, `burnin`, `thin`. Sampling is uniform over
all eligible records, with replacement by default. See [posterior input](trees.md#posterior-tree-streams-and-samples).

### Posterior distance means

`mean_dist(data, trees; index_family=:sorensen, component=:total, pairs=nothing,
nan_policy=:skip, block_size=256)` consumes trees directly. Matrix columns
without `tip_labels` take their meaning from the first tree. `CommunityData`
supplies species labels. The one-argument form aggregates existing distances
or beta results, aligning site labels.

`MeanDistResult` contains `sample_labels`, `distances`, `counts`, `pairs`,
`replicates`, `component`, `index_family`. `counts` is per pair and counts only
finite observations; `replicates` counts all inputs. No valid value gives NaN.
`nan_policy=:propagate` preserves NaN if any replicate is undefined. See
[posterior means](regionalization.md#mean-distances-across-posterior-trees).

### Grouping and conservation results

| Function | Result and ordering | Guide |
| --- | --- | --- |
| `optimal_functional_groups` | `FunctionalGroupSelectionResult`: sorted candidate k, diagnostics, optional recommendation and reusable grouping | [Aggregation](aggregation.md) |
| `nmds` | `NMDSResult`: site-order `points`, Kruskal `stress`, packed input/score/fitted distances, best-start and per-start diagnostics | [NMDS](regionalization.md#nonmetric-multidimensional-scaling) |
| `align_nmds` | `NMDSAlignmentResult`: source-order coordinates, matched labels, rotation/scale/centers, anchor residuals and RMSE | [Alignment](regionalization.md#coordinate-alignment-and-colors) |
| `nmds_colors` | `NMDSColorResult`: labels, two-axis `points`, `rgb`/`hex`, reusable `NMDSColorScale`, `clipped` | [Colors](regionalization.md#coordinate-alignment-and-colors) |
| `fitgom` | `GoMResult`: labeled sites × k `omega`, species × k `theta`, hard memberships, objective and convergence diagnostics | [Mixed membership](regionalization.md#mixed-regional-membership) |
| `indicators` | `IndicatorSpeciesResult`: regions × ranks, original species indices, missing padding | [Indicators](regionalization.md#regional-indicator-species) |
| `phylobuilder` | `PhyloBuilderResult`: independent tree; status and original anchor IDs in requested-species order | [Trees](trees.md) |
| `edge` | `EDGEResult`: full-tree tip order, ED/GE/EDGE and missing-category status | [Conservation](conservation.md) |
| `hotspots`, `coldspots` | `MetricSelectionResult`: input-order flags and labels, quantile threshold and valid counts | [Conservation](conservation.md) |
| `classify_values` | `ValueClassificationResult`: input-order classes, upper bounds and optional Jenks SSE | [Conservation](conservation.md) |

GoM weights describe composition, not uncertainty. Convergence of GoM or NMDS
does not establish a global optimum. Metric-value summaries allow finite
negative numbers; `invalid=:skip` keeps missing/NaN positions as missing
flags/classes. Community abundances still must be nonnegative.

### Trait Result Types

| Function | Result and main fields |
| --- | --- |
| `pagel_lambda` | `PagelLambdaResult`: fitted lambda, profile interval, likelihood-ratio test, fit diagnostics and model summary; DataFrame input aligns taxon labels |
| `phylogenetic_signal_k` | `PhyloSignalResult`: `k`, observed/null PIC variance, null SD, rank p-value, z score, aligned traits/tip labels, valid count and status |
| `phylogenetic_contrasts` | `PICResult`: contrasts, internal-node IDs, expected distances, local ancestral values, tip labels, `scaled`/`weighted` |

K shuffles tip traits; its default RNG is `StableRNG(20260713)`. PIC requires
a strictly bifurcating tree with finite positive branches; use `multi2di`
explicitly for multifurcations. PIC is not a general ancestral-state
reconstruction API. See [trait models and diagnostics](signals.md).

### Boundary Policies

| Case | Behavior |
| --- | --- |
| Empty site / singleton PD | `pd`: zero for empty sites; root-to-tip length for singletons with `include_root=true`, NaN with false. `pd_calc` variants have their own branch rules; singleton MST is zero. |
| Empty MPD/MNTD | NaN; singletons also NaN except abundance-weighted MPD with positive mass, which is zero |
| Fewer than two species | PSV/PSR/PSE/PSC and PSV/PSR variances are NaN |
| Empty/zero-PD beta | `beta_diss` uses richness, `phylobeta` uses PD. Both zero: all components NaN. Only one zero: total one, turnover/nestedness NaN. |
| Zero UniFrac/PhyloSor denominator | Both return zero, including nonempty sites whose PD is zero |
| Invalid community values | Negative, infinite or NaN values raise `ArgumentError`; rows are never silently dropped |
| Tree restrictions | RaoD requires ultrametricity (absolute root-to-tip tolerance `1e-8`); scaled PSV-family metrics require positive root-to-tip lengths |
| Missing/nonfinite traits | K/PIC reject them |
| Randomized methods | Pass an explicit `AbstractRNG` for reproducibility. `nri_nti` also accepts stable seeds; defaults differ across functions. |

SES zero-SD behavior and lower-tail p-value denominators are specified in
[null models](null_models.md). Clustering/NMDS require complete, finite,
nonnegative symmetric distances with zero diagonal; see [regionalization](regionalization.md).
PhyloSor is a similarity, comdist diagonals are generally nonzero, and SES
z scores are not dissimilarity inputs.
