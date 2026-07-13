# EvoKits Public API Contract

This document freezes the public conventions for the 0.2 API. All functions
accept a `CompactTree` unless the signature below says otherwise.

## Shared Ordering

`tree.tip_labels` is the canonical tip order. A tip position is an integer in
`1:tree.ntips`; a compact node ID is an integer in `1:tree.nnodes`.
Functions that accept labels resolve them through `tree.tipname_to_id`.

Community matrices use rows as communities and columns as tips. The row order
is preserved in every returned vector or matrix. `tip_labels` describes the
input column order and reorders columns into `tree.tip_labels` order; it does
not rename or reorder community rows. Duplicate, missing, or unknown labels
raise `ArgumentError`.

Community values must be finite and non-negative. Values greater than zero
indicate presence; abundance-weighted methods use the positive values as
weights.

Trait vectors are aligned to tree order by default. When `tip_labels` is
provided, it describes the order of the supplied trait vector and the vector
is reordered to `tree.tip_labels`. Missing and non-finite trait values are
rejected by K and PIC.

## Matrix Semantics

`comdist`, `comdistnt`, `unifrac`, and `phylosor` return ordinary `n x n`
`Float64` matrices. Matrix indices follow community row order.

- `comdist`: mean pairwise distance across the two communities. The diagonal
  is within-community mean distance; an empty community gives `NaN`.
- `comdistnt`: symmetric mean of the two directional nearest-taxon means. The
  diagonal is the within-community nearest-taxon mean. With
  `exclude_conspecifics=true`, a singleton diagonal is `NaN`.
- `unifrac`: unweighted branch-length distance. The diagonal is zero,
  including an empty-versus-empty comparison; an empty-versus-nonempty pair
  is one when the nonempty community spans positive branch length.
- `phylosor`: `2 * shared_PD / (PD1 + PD2)`. The diagonal is one for a
  nonempty community and zero for an empty community.

`raoD` returns `RaoDResult` with `Dkk`, `Dkl`, `H`, `total`, `alpha`, `beta`,
and `Fst`, matching Picante's partitioning convention. It requires an
ultrametric tree and uses distance-to-MRCA (`patristic / 2`).

## Community Result Types

- `pd` returns `PDResult` with `sample_labels`, `pd`, `richness`, and
  `include_root`.
- `ses_pd`, `ses_mpd`, and `ses_mntd` return typed results containing observed
  values, null means/SDs, z scores, rank p-values, run count, and null model.
  Empty or singleton MPD/MNTD observations are `NaN`.
- `psv` returns `PSVResult` with `psv`, `richness`, `variance`, and
  `scale_vcv`. Communities with fewer than two taxa have `NaN` PSV and
  variance. `psr`, `pse`, `psc`, and `psd` return named tuples whose field names
  are the metric names documented by their function.
- `expected_pd` and `pd_variance` return `ExpectedPDResult`; `richness` is
  `1:tree.ntips`, `expected_pd` is the analytical expectation, and
  `variance_pd` is the upper-bound or exact Picante variance selected by
  `upper_bound`.
- `pd_calc(tree, tips; variant=...)` returns one scalar for a tip subset.
  Matrix input returns a named tuple with `pd`, `richness`, and normalized
  `variant`. Accepted variants are `TBL`, `MST`, `UEH`, `SBL`, and `TIP`.
- `pd_bootstrap(tree, ntips; ...)` returns a vector of replicate scalars.
  `pd_bootstrap(community, tree; ...)` returns `PDBootstrapResult` with
  `observed`, `replicates`, and `richness`.

## Trait Result Types

`phylogenetic_signal_k` returns `PhyloSignalResult` with `k`, observed and
null PIC variance, null SD, rank p-value, z score, reordered trait values,
tip labels, valid replicate count, and status. Its null model is tip shuffling;
pass an explicit `rng` for reproducibility. The default is `StableRNG(20260713)`.

`phylogenetic_contrasts` returns `PICResult` with `contrasts`, internal-node
IDs, expected contrast distances, local ancestral values, tip labels, and the
`scaled`/`weighted` flags. PIC requires a strictly bifurcating tree with finite
strictly positive branch lengths. Use `multi2di` explicitly for multifurcating
trees. PIC is not a general ancestral-state reconstruction API.

## Boundary Policies

- `pd` assigns zero PD to an empty community and the root-to-tip path length
  to a singleton when `include_root=true`.
- Empty or undefined community metrics return `NaN` rather than silently
  dropping rows; UniFrac and PhyloSor use the explicit diagonal rules above.
- Negative, infinite, or NaN community values raise `ArgumentError`.
- K/PIC reject missing and non-finite trait values.
- RaoD rejects non-ultrametric trees. PIC rejects multifurcations and zero or
  non-finite branch lengths.
- Randomized methods accept an injected `AbstractRNG`; seed-based
  `nri_nti` remains deterministic for the same seed.

The package deliberately does not export BM/OU/Mk fitting, general ASR,
PGLS, diversification, dating, SIMMAP, plotting, or Castor C++ wrappers.
