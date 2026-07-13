# EvoKits Structural and Algorithm Migration Plan

## 1. Objective

EvoKits will have exactly four source subdirectories:

    src/
      EvoKits.jl
      base/
      tree/
      traits/
      community/

There will be no src/tools/, src/core/, or nested src/tree/mapping/
directory after the restructuring.

The algorithm scope is fixed to all P0 algorithms, all P1 algorithms, and
only the two P2 methods K signal and Felsenstein PIC.

The restructuring must preserve the existing behavior of CompactTree, PD/SES,
NRI/NTI, Pagel lambda, EDGE2, Phylo.D, pruning, and PhyloRef before any new
algorithm is added.

## Current Implementation Checkpoint

The four-layer refactor and P0 tree kernel are complete. `CompactTree` and its
Newick boundary now live in `base/tree_types.jl` and `base/tree_io.jl`;
tree-only algorithms remain in `tree/`.

P1 currently includes the shared community input contract, beta diversity
(`comdist`, `comdistnt`, UniFrac, PhyloSor, RaoD), PSV-family metrics, edge
abundance distribution, expected PD/variance, explicit `pd_calc` variants,
and tip-weighted PD bootstrap. P2-restricted K signal and Felsenstein PIC are
implemented with injected RNG, tip-label alignment, and explicit branch-length
and multifurcation validation. Remaining work is reference-fixture expansion,
large-tree benchmarks, and final API review; the current Picante/Caper
hand-fixture and boundary suite is recorded in `test/runtests.jl`, with scale
results in `PERFORMANCE_BENCHMARKS.md`. Excluded model layers remain out of
scope.

The current verification checkpoint is complete for this migration pass:

- P1 Picante/Caper hand fixtures and edge-case policies pass.
- P2 K/PIC hand calculations and invalid-input policies pass.
- 1k/10k/100k-tip smoke benchmarks are recorded.
- Public exports contain no excluded model or workflow API.

## 2. Scope Lock

### Allowed

- Tree topology, I/O, identity, distances, MRCA, clades, pruning, and tree
  shape statistics.
- Community PD, MPD, MNTD, SES, beta diversity, UniFrac, PhyloSor, RaoD,
  PSV-family metrics, expected PD, PD variants, and ED corrections.
- Local-pool NRI/NTI, EDGE2, and existing conservation metrics.
- Continuous-trait K signal and Felsenstein PIC.

### Prohibited

Do not implement, re-export, or create placeholder APIs for:

- BM, OU, ROU, Mk, multi-Mk, threshold, hidden-state, or stochastic-mapping
  models.
- General ASR, fastAnc, Bayesian ASR, maximum-parsimony ASR, or squared-change
  ASR. PIC may calculate only the ancestral values needed by its contrasts.
- PGLS, CAIC, crunch, brunch, macrocaic, or comparative.data.
- Diversification, HBD/HBDS, MuSSE/DSSE, gene-tree, dating, fossil, or
  past-diversification workflows.
- Plotting, graphics, SIMMAP visualization, R compatibility objects, or Castor
  C++ wrappers.

No excluded feature should receive a placeholder export or an early partial API.

## 3. Final Module Boundaries

### 3.1 base

Package-wide foundations and the shared tree data boundary:

    base/criteria.jl      AIC, AICc, BIC
    base/random.jl        stable seed and RNG conventions
    base/tree_types.jl    CompactTree storage representation
    base/tree_io.jl       Newick <-> CompactTree conversion and file I/O

`CompactTree` and its Newick conversion boundary are intentionally placed in
base because every higher layer consumes the same immutable storage contract.
Tree algorithms, CompactTree validation, topology queries, pruning, and shape
statistics remain in tree. Community matrix alignment and null models remain
domain-specific and do not belong in base.

### 3.2 tree

All tree representation, topology, identity, and tree-only statistics:

    tree/validation.jl
    tree/query.jl
    tree/distances.jl
    tree/clades.jl
    tree/pruning.jl
    tree/mapping.jl
    tree/transforms.jl
    tree/multi2di.jl
    tree/metrics.jl
    tree/shape.jl

Keep types, io, and pruning in tree. Merge tree/mapping.jl and its nested
files into the single tree/mapping.jl. Move the LCA index, virtual-tree
queries, and tree path helpers out of tools/community_metrics.jl and
tools/common.jl into query.jl or clades.jl. Move make_ultrametric.jl and
multi2di.jl into transforms.jl. Move evol_distinct.jl into metrics.jl.

The tree layer must not depend on traits or community.

### 3.3 traits

Trait alignment, existing comparative utility, and the two allowed P2 methods:

    traits/alignment.jl
    traits/pagel_lambda.jl
    traits/phylo_d.jl
    traits/signal_k.jl
    traits/pic.jl

Move the current alignment file here. Move pagel_lambda.jl and phylo_d.jl
from tools here. Phylo.D remains a binary-trait signal method, not a general
discrete-trait model.

Traits may depend on base and tree, but not community.

### 3.4 community

Community input, null models, community metrics, and conservation metrics:

    community/input.jl
    community/null_models.jl
    community/pd.jl
    community/mpd.jl
    community/mntd.jl
    community/beta.jl
    community/diversity.jl
    community/nri_nti.jl
    community/edge2.jl

Split community input and randomization code from tools/common.jl. Move
ses_pd.jl, ses_mpd.jl, ses_mntd.jl, nri_nti.jl, and edge2.jl here. Extract
the tree-query portion of tools/community_metrics.jl before deleting it.

Community may depend on base and tree, but tree must never depend on community.

## 4. Include Order

EvoKits.jl should contain dependencies, exports, and includes only. Use this
order:

    base/criteria.jl
    base/random.jl
    base/tree_types.jl
    base/tree_io.jl
    tree/validation.jl
    tree/query.jl
    tree/distances.jl
    tree/clades.jl
    tree/pruning.jl
    tree/mapping.jl
    tree/transforms.jl
    tree/multi2di.jl
    tree/metrics.jl
    tree/shape.jl
    traits/alignment.jl
    traits/pagel_lambda.jl
    traits/phylo_d.jl
    traits/signal_k.jl
    traits/pic.jl
    community/input.jl
    community/null_models.jl
    community/pd.jl
    community/mpd.jl
    community/mntd.jl
    community/beta.jl
    community/diversity.jl
    community/nri_nti.jl
    community/edge2.jl

The order must not create circular dependencies. Do not add a directory-level
dispatcher.

## 5. Phase 0: Structural Refactor Without New Algorithms

1. Snapshot the current public exports and complete test results.
2. Move files with history-preserving git operations where possible.
3. Split tools/common.jl into base, tree, and community responsibilities.
4. Extract tree query code from tools/community_metrics.jl.
5. Flatten tree/mapping/ into sibling mapping files in tree/ (or merge them
   into mapping.jl once the history-preserving move is complete).
6. Move CompactTree-specific validation into tree/validation.jl.
7. Rename files into the four directories and rewrite EvoKits.jl includes.
8. Preserve all public names and numerical behavior.

Gate:

- src contains only EvoKits.jl plus base, tree, traits, and community.
- using EvoKits succeeds in a clean environment.
- Existing tests pass with unchanged results.
- tree has no dependency on community or traits.
- git diff --check passes.

No new P0/P1 algorithm starts before this gate passes.

## 6. Phase P0: Tree Kernel and Tree Shape

### 6.1 Query and distance kernel

Add to tree/query.jl, distances.jl, and clades.jl:

    mrca(tree, a, b)
    mrca(tree, tips)
    distance(tree, a, b; weighted = true)
    pairwise_distance(tree; tips = nothing, weighted = true)
    distances_from_root(tree; weighted = true)
    distances_to_tip(tree, focal; weighted = true)
    descendant_tips(tree, node; labels = false)
    clade_members(tree; include_nodes = false)
    subtree_at_node(tree, node)

Freeze accepted node forms: compact node IDs, tip positions, and labels.
Document matrix order and diagonal behavior. Reuse the binary-lifting index,
child-to-edge index, virtual-tree routines, and shared-root covariance logic.

### 6.2 Tree shape

Add to tree/shape.jl:

- Sackin index.
- Colless and normalized Colless indices.
- Blum index.
- Pybus gamma for rooted ultrametric trees.
- Lineage-through-time counts.
- Node depth, node height, and branch/path summaries.

Caper Fusco imbalance is considered only after confirming it is distinct from
the selected imbalance statistics.

### 6.3 P0 tests

- Balanced, ladder, star, non-ultrametric, and multifurcating fixtures.
- Brute-force path comparison for every pairwise distance.
- MRCA symmetry and ancestor identity.
- Weighted distances equal edge-length sums; unweighted distances equal edge
  counts.
- Subtree tip order and branch-length preservation.
- Known imbalance values and gamma error cases.
- 1k, 10k, and 100k-tip performance measurements.

## 7. Phase P1: Community Algorithms

### 7.1 Shared input and null-model contract

Centralize in community/input.jl:

- Rows are communities and columns are tips.
- Positive values indicate presence.
- Abundance-weighted methods use positive values as weights.
- tip_labels reorders columns to tree-tip order.
- Empty and singleton behavior is explicit.

Centralize in community/null_models.jl:

- Taxa-label, richness, frequency, sample-pool, phylogeny-pool,
  independent-swap, and trial-swap randomization.
- Injected RNG, stable rank p-values, online null statistics, and explicit
  zero-variance status.

Existing PD, MPD, MNTD, SES, and NRI/NTI behavior must remain unchanged.

### 7.2 Beta diversity

Implement in community/beta.jl:

    comdist(community, tree; abundance_weighted = false, tip_labels = nothing)
    comdistnt(community, tree; abundance_weighted = false,
              exclude_conspecifics = false, tip_labels = nothing)
    unifrac(community, tree; tip_labels = nothing)
    phylosor(community, tree; tip_labels = nothing)
    raoD(community, tree; abundance_weighted = true, tip_labels = nothing)

Use the P0 distance API rather than calculating paths inside each metric.
Return labeled symmetric matrices with explicit diagonal semantics.

Test identical, disjoint, nested, singleton, empty, weighted, and unweighted
communities against hand calculations and Picante fixtures.

### 7.3 Diversity and PD extensions

Implement in community/diversity.jl:

- PSV, PSR, PSE, PSC, and PSD.
- Edge-abundance distribution.
- Analytical expected PD and PD variance.
- Caper pd.calc variants TBL, MST, UEH, SBL, and TIP.
- PD bootstrap with optional tip weights.
- Isaac, Mooers, and no-correction ED polytomy modes.

Do not silently change the current pd default. Add explicit keywords or clearly
named new functions for alternative semantics. Define root-edge, singleton,
empty-community, zero-length, and singular-covariance policies.

### 7.4 P1 tests

- Enumerate all small-tree subsets to verify expected PD and variance.
- Hand-calculate beta-diversity examples.
- Symmetry and range properties for UniFrac and PhyloSor.
- Reference comparisons for all PD variants and ED corrections.
- Null-model reproducibility and degenerate-null tests.
- Allocation benchmarks for large community matrices.

## 8. Phase P2-Restricted: K Signal and PIC

### 8.1 K signal

Add traits/signal_k.jl:

    PhyloSignalResult
    phylogenetic_signal_k(tree, trait;
                          permutations = 999,
                          rng = StableRNG(20260713),
                          tip_labels = nothing)

Implement only the Picante/Phytools K definition and tip-shuffle test. Return
K, observed PIC variance, null mean and SD, rank p-value, z-score, trait
labels, and valid replicate count. Reuse shared covariance and factorization
helpers. Do not add lambda or any other model fitting.

### 8.2 PIC

Add traits/pic.jl:

    PICResult
    phylogenetic_contrasts(tree, trait;
                           scaled = true,
                           weighted = true,
                           tip_labels = nothing)

This is Felsenstein independent contrasts only. It is not a general ASR
engine. Freeze multifurcation, zero-length, missing-value, and non-ultrametric
semantics before implementation. Initially reject missing values and document
multi2di as the route for analyses requiring a binary tree.

Test hand-calculated two-tip and four-tip contrasts, reference fixtures,
scaled/unscaled modes, non-ultrametric trees, invalid branch lengths, missing
values, and misaligned labels.

## 9. Documentation and API Cleanup

After each phase:

- Update README examples and exports.
- Document dimensions, label order, RNG, null model, and missing-value policy.
- Remove internal _tools_ names from user-facing documentation.
- Keep compatibility aliases only when they do not duplicate algorithms.
- Add one end-to-end example for P0, P1, K, and PIC.
- Describe the four source layers rather than the historical extraction layout.

## 10. Verification Gates

Every phase must pass:

1. using EvoKits in a clean environment.
2. Pkg.test() with unchanged existing numerical results.
3. Reference fixtures and independent brute-force tests.
4. git diff --check.
5. A recorded performance benchmark at the intended scale.
6. No inverted dependency between the four layers.
7. No excluded algorithm added indirectly as a helper or compatibility wrapper.

## 11. Execution Order

1. Baseline tests and public API snapshot.
2. Four-directory structural refactor with no semantic changes.
3. P0 query, distance, clade, and subtree kernel.
4. P0 shape, gamma, and LTT statistics.
5. P1 community input and null-model consolidation.
6. P1 beta diversity.
7. P1 PSV, expected PD, PD variants, bootstrap, and ED corrections.
8. P2 K signal.
9. P2 PIC.
10. Documentation, performance reports, and final export review.

No work outside this list should start without explicitly revising this plan.
