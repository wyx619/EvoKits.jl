# API conventions

[Home](../README.md) · [简体中文](api_zh.md)

A combined example, ordering and result contracts, and boundary behavior.

## End-to-End Example

This example combines tree, community and trait analysis in one session:

```julia
using EvoKits
using Random

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end

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

### Posterior tree input

`read_trees(paths; burnin=0, thin=1)` returns a single-pass `TreeStream` of
`CompactTree`s; its do-block closes files on early exit. `tree_sampler` returns
an iterable `TreeSampleResult`, with independent trees and draw-order indices,
source paths, per-file record numbers and line numbers. Filtering is per file;
sampling is uniform across all eligible records, with replacement by default.
Newick/NEXUS detection is per file. NEXUS indices count only `TREE`/`UTREE`
definitions; translation tables apply within their own `TREES` block.
See [tree input](trees.md#posterior-tree-streams-and-samples) for format and error rules.

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

`community_longtable(data)` exports a named tuple of `table`, `sample_labels`,
and `species_labels`. The table has one row per positive value, with `sample`,
`species`, and `value` columns in species-then-site order. Preserve both full
label lists when reconstructing with `CommunityData(table; sample_labels,
species_labels)` to retain empty domains and original order. Values preserve
abundance and numeric type; reconstruction follows long-table rules, including
Bool-to-Int counts. See [the round-trip example](community.md#export-a-long-table-and-reconstruct-the-community).

`comdist`, `comdistnt`, `unifrac`, and `phylosor` return ordinary `n x n`
`Float64` matrices by default. Matrix indices follow community row order.

- `comdist`: mean pairwise distance across the two communities. The diagonal
  includes zero self-distances in the within-community mean; an empty site gives `NaN`.
- `comdistnt`: equal mean of the two directional nearest-taxon means, or a mean
  weighted by original site totals with `abundance_weighted=true`. Nonempty
  diagonals are zero by default. With `exclude_conspecifics=true`, diagonals
  are within-site nearest-other-species means (a singleton gives `NaN`).
  Exclusion uses species identity; distinct species at distance zero remain eligible.
- `unifrac`: unweighted branch-length distance. The diagonal is zero,
  including an empty-versus-empty comparison; an empty-versus-nonempty pair is
  one when the nonempty community spans positive branch length.
- `phylosor`: `2 * shared_PD / (PD1 + PD2)`. The diagonal is one for a
  nonempty community and zero for an empty community.

Both `comdist` and `comdistnt` accept `pairs`, returning a `Vector{Float64}`
in request order with repeated, reversed and self pairs retained. Indices
refer to the complete site order; an empty request returns an empty vector.
All input values, labels and branch lengths are validated.

`prepare_community_distances(matrix_or_data, tree; pairs=nothing)` creates a
`CommunityDistanceContext` for both distance functions and both weighting
modes. Matrix input also accepts `tip_labels` and `sample_labels`; site labels
must be unique and nonempty. Without preparation-time pairs, context calls
default to full matrices and allow changing `pairs`. Preparation-time pairs
fix the vector request for later calls, retaining only its effective endpoints;
changing that request requires rebuilding. Input edits never update the snapshot.
Returned arrays are independent; shared-context calls are serialized.
See [the repeated-call example](community.md#distances-between-communities).

`raoD` returns `RaoDResult` with `Dkk`, `Dkl`, `H`, `total`, `alpha`, `beta`,
and `Fst`. It uses `patristic / 2` and requires an ultrametric tree (absolute
root-to-tip tolerance 1e-8), with finite non-negative branch lengths and a
finite sum. Abundance weighting defaults to true. Setting it false changes
within-site compositions to incidence; pooling and site weights retain original
abundances. Empty sites propagate NaN to alpha/beta/Fst while pooled total may
remain valid. Dkk, Dkl and H own independent arrays in input site order.
See [Rao partition rules and example](community.md#rao-diversity-partition).

### Community Result Types

- `pd` returns `PDResult` with `sample_labels`, `pd`, `richness`, and
  `include_root`.
- `weighted_endemism` returns `WeightedEndemismResult` with `sample_labels`,
  `we`, `richness`, and `species_occupancy` in input column order.
- `phylo_endemism` returns `PhyloEndemismResult` with `sample_labels`, `pe`,
  `richness`, `branch_occupancy` in tree edge order, and `weighted`.
- `relative_pd` returns `RelativePDResult` with `sample_labels`, `rpd`,
  `pd_normalized`, `pd_equal_normalized`, `richness`, `total_branch_length`,
  and `positive_branch_count`.
- `phylobeta` returns `PhylobetaResult` with `sample_labels`, `total`,
  `turnover`, `nestedness`, `index_family`, and `pairs`.
- `beta_multi` returns `BetaMultiResult` with scalar `total`/`turnover`/`nestedness`,
  `index_family`, selected `sample_labels`, original `sample_indices`, selected-site
  `richness` and `pooled_richness`. See [multiple-site beta](community.md#multiple-site-taxonomic-beta).
- `phylobeta_multi` returns `PhylobetaMultiResult` with the three scalar components,
  `index_family`, selected `sample_labels`, original `sample_indices`, `richness`,
  rooted `pd` and selected-union `pooled_pd`. It also accepts `CommunityContext`.
  See [multiple-site phylogenetic beta](community.md#multiple-site-phylogenetic-beta).
- `phylobeta_ses` returns `PhylobetaSESResult` with site labels and richness,
  `observed`, `null_mean`, `null_sd`, `z_score`, `p_value`,
  `valid_replicates`, `pairs`, the selected family/component, and null settings.
- `nmds` returns `NMDSResult` with labels, `points`, Kruskal `stress`, packed
  input/score/fitted distances, the best start's history and stopping reason,
  and per-start stress, iteration and convergence diagnostics.
- `align_nmds` returns `NMDSAlignmentResult` with source-order coordinates,
  matched labels, the rotation/scale/centers, anchor residuals and RMSE.
- `nmds_colors` returns `NMDSColorResult` with labels, two-axis `points`,
  matching `rgb`/`hex`, a reusable `NMDSColorScale` and `clipped` flags.
  See [alignment and colors](regionalization.md#coordinate-alignment-and-colors).
- `fitgom` returns `GoMResult` with labeled `omega`/`theta`, hard membership,
  the NEF objective, EM residual and per-start convergence diagnostics.
  `indicators(fit)` ranks its species. See [mixed membership](regionalization.md#mixed-regional-membership).
- `ses_pd`, `ses_mpd`, and `ses_mntd` return typed results containing observed
  values, null means/SDs, z scores, rank p-values, run count, and null model.
  Empty communities have `NaN` MPD/MNTD. Singleton observations are `NaN`,
  except abundance-weighted MPD with positive mass, which is zero.
- `psv` returns `PSVResult` with `psv`, `richness`, `variance`, and `scale_vcv`.
  Variance describes uniform sampling without replacement at that richness
  from all tree tips; full-pool variance is zero. This corrects the previous
  leaf-order-dependent formula. Communities with fewer than two taxa have
  `NaN` PSV and variance; `compute_var=false` leaves variance as `NaN`.
  `psr`, `pse`, `psc`, and `psd` return named tuples whose field names are
  the metric names documented by their function. See [PSV and PD moments](community.md#psv-psr-and-analytical-pd-moments).
- `pse` uses relative abundances; equal positive abundances give PSV. `psc`
  uses presence and excludes the same species identity, retaining zero
  contrasts to different species. Both return their metric and `richness`,
  with `NaN` for fewer than two species. Raw-covariance PSC can be negative.
  `psd` returns `psv`, `var_psv`, `psc`, `psr`, `var_psr`, `pse` and `richness`
  in site order, each with its own array and the individual metric's rules.
- `expected_pd` and `pd_variance` return `ExpectedPDResult`; `richness` is
  `1:tree.ntips`. Every tip independently occurs with probability
  `richness/tree.ntips`, so richness is an expected count. `expected_pd` is
  the mean PD; `variance_pd` is the exact Bernoulli variance when
  `upper_bound=false`, or picante's covariance upper bound by default.
- `pd_calc(tree, tips; variant=...)` returns one scalar for a tip subset.
  Matrix input returns a named tuple with `pd`, `richness`, and normalized
  `variant`. Accepted variants are `TBL`, `MST`, `UEH`, `SBL`, and `TIP`.
- `pd_bootstrap(tree, ntips; ...)` returns a vector of replicate scalars.
  `pd_bootstrap(community, tree; ...)` returns `PDBootstrapResult` with
  `observed`, `replicates`, and `richness`.

### Grouping and conservation results

| Function | Result and ordering | Guide |
|---|---|---|
| `optimal_functional_groups` | `FunctionalGroupSelectionResult`: sorted candidate k, diagnostics, optional recommendation and reusable grouping | [Aggregation](aggregation.md) |
| `indicators` | `IndicatorSpeciesResult`: regions × ranks, original species indices, missing padding | [Regionalization](regionalization.md) |
| `phylobuilder` | `PhyloBuilderResult`: independent tree; status and original anchor IDs in requested-species order | [Trees](trees.md) |
| `edge` | `EDGEResult`: full-tree tip order, ED/GE/EDGE and missing-category status | [Conservation](conservation.md) |
| `hotspots`, `coldspots` | `MetricSelectionResult`: input-order flags and labels, quantile threshold and valid counts | [Conservation](conservation.md) |
| `classify_values` | `ValueClassificationResult`: input-order classes, upper bounds and optional Jenks SSE | [Conservation](conservation.md) |

These operations preserve their inputs. Vector metric summaries accept finite
negative values; `invalid=:skip` preserves missing/NaN positions with missing
flags/classes. This policy is separate from nonnegative community inputs.

### Trait Result Types

`pagel_lambda` returns `PagelLambdaResult` with the fitted lambda, profile
interval, likelihood-ratio test, fit diagnostics, and model summary fields
described in the [phylogenetic signal guide](signals.md). The DataFrame
overload performs taxon-label alignment before fitting.

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
