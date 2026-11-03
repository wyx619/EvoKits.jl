# Regionalization and ordination

[Home](../README.md) · [简体中文](regionalization_zh.md)

Cluster beta distances, fit mixed memberships, evaluate partitions, average posterior distances, run NMDS and rank regional indicator species.

## Regionalization from beta diversity

`phyloregion` clusters a complete distance matrix and cuts the full hierarchy
into an explicit number of regions `k`, using UPGMA by default. It also
accepts full `BetaDissResult` and `PhylobetaResult` objects, retaining site labels.
Select `component=:total` (default), `:turnover`, or `:nestedness`.

```julia
using EvoKits, SparseArrays

regional_data = CommunityData(sparse([1 1 0 0; 1 1 1 0; 0 0 1 1; 0 1 1 1]);
    species_labels = ["A", "B", "C", "D"],
    sample_labels = ["west_1", "west_2", "east_1", "east_2"])
beta = beta_diss(regional_data)
regions = phyloregion(beta; k = 2)
regions.membership       # input-row order; region IDs start at 1
regions.cluster_sizes
regions.region_distances
regions.distinctiveness

# Or supply a complete distance matrix and its labels.
regions = phyloregion([0 1 8 10; 1 0 6 12; 8 6 0 2; 10 12 2 0];
    k = 2, sample_labels = ["A", "B", "C", "D"])
# membership = [1, 1, 2, 2]; region_distances = [0 9; 9 0]
```

`PhyloregionResult` returns site `sample_labels` and `membership`,
`cluster_sizes`, `region_distances`, `distinctiveness`, and the complete
hierarchy (`merge` and `heights`), along with `k`, `method`, and source
`component`. Regions are numbered by their first input site. Negative merge
entries identify input rows; positive entries identify earlier merge rows.
Region distances average all cross-site pairs between two regions and have
zero diagonals. Distinctiveness averages distances to the other regions with
equal region weight, regardless of region size; it is NaN for k=1.

Input must be nonempty, square, exactly symmetric, finite, non-negative, with
zero diagonal; it need not be Euclidean. Float64 conversion must preserve
finite positive distances. Undefined values, including NaN diagonals from
empty sites, are rejected; explicitly select a defined subset and its labels
before calling. No sites are automatically removed or distances imputed.
Selected-pair beta results cannot supply the complete distances needed here.

Supported methods follow R `hclust` recurrences:

| `method` | Meaning |
|---|---|
| `:average` | UPGMA; average of all cross-site distances |
| `:single` | Minimum cross-site distance |
| `:complete` | Maximum cross-site distance |
| `:mcquitty` | WPGMA; equal weight for the two child clusters |
| `:ward_D` | R's older Ward recurrence on supplied distances |
| `:ward_D2` | Ward recurrence on squared distances, square-root heights |
| `:centroid` | UPGMC recurrence on supplied distances |
| `:median` | WPGMC recurrence on supplied distances |

`Symbol("ward.D")` and `Symbol("ward.D2")` are accepted aliases; results use
`:ward_D` and `:ward_D2`. Ward's variance interpretation requires Euclidean
distances. Centroid/median do not automatically square the input; their heights
can decrease and remain in chronological merge order. A linkage update outside
finite Float64 range is an error; rescale input distances in that case.

The six monotone methods use O(n²) nearest-neighbor chaining. Centroid/median
use global minima with cached row minima, with O(n³) worst-case time. Both use
a packed Float64 work triangle of `n*(n-1)/2` entries. The full input and
regional output still require quadratic memory. Ties are deterministic for
fixed row order: chains start at the earliest active cluster, retain the
predecessor on equal nearest distances, otherwise choose the earliest active
cluster; equal-height discoveries keep their order. Centroid/median choose the
lexicographically first active pair. Ties can differ from R `hclust` and depend
on row order. Inputs are not mutated; outputs own their arrays.

### Reuse a hierarchy and evaluate cuts

```julia
ctx = prepare_phyloregion(beta_diss(regional_data); method = :average)
cutree(ctx; k = 2)                 # memberships only; O(n), no reclustering
regions3 = phyloregion(ctx; k = 3) # regional summaries, reusing the same tree
cutree(regions3; k = 2)            # any PhyloregionResult retains the full tree

selection = select_linkage(ctx)    # compare all eight methods
selection.methods
selection.correlations            # original vs cophenetic Pearson correlation
selection.best_method             # first finite maximum, or nothing

evaluation = optimal_phyloregion(ctx; k = 4) # evaluate every cut from 1 to 4
evaluation.df                     # columns: k, ev, totwss, totbss, tss
evaluation.optimal_k               # suggested k, or nothing
evaluation.reason
```

`PhyloregionContext` owns an original packed distance snapshot and full tree;
input changes do not affect it. Treat context arrays as read-only. Preparation
retains one packed triangle and uses a second as temporary clustering
workspace. Repeated `cutree` calls only cut the tree; `phyloregion(ctx; k)` also
visits site pairs for regional means and allocates the k×k regional matrix.

`select_linkage` also accepts a distance matrix or full beta result. Restrict
methods with, for example, `methods=(:average, :complete, :single)`. It compares
original distances with the height where each pair first joins the tree.
Larger correlation means better distance preservation, not ecological
correctness. Constant distances/heights and a single pair give NaN; if all
scores are undefined, `best_method` is `nothing`. No full cophenetic matrix is
allocated. With a context, its existing tree is reused for the matching method.

`optimal_phyloregion` also accepts a distance matrix or full beta result,
building one hierarchy. `k` is the upper bound, defaulting to `min(n,20)`.
Within-cluster sums are the sum of squared unordered distances divided by the
cluster's site count; TSS uses all sites, BSS=TSS−WSS, and `ev=BSS/TSS`.
For non-Euclidean dissimilarities these are distance-based scores; EV may be
nonmonotone or outside [0,1].

Selection follows the saved phyloregion threshold heuristic: try
`ev_thresholds=(0.95,0.9,0.8,0.75,0.67,0.5,0.33)` descending, then
`inc_thresholds=(0.01,0.05,0.1)` ascending. Find the first k meeting each EV
threshold and the first whose next EV increment is small enough, and recommend
their maximum for the first successful threshold combination. Tolerance is
`10^(-precision)`, default precision=3; the last evaluated k has no next
increment. This is a heuristic, not a curvature optimization. Unmet thresholds
or zero dispersion return `optimal_k=nothing` and the full table, with
`:thresholds_not_met` or `:zero_dispersion` as the reason. Sums use original
squared-distance units and can overflow/underflow, while EV and decisions use
normalized distances. Scores visit all pairs once with O(n) extra storage.

For ordination coordinates, see [NMDS](#nonmetric-multidimensional-scaling).

## Mixed regional membership

`fitgom` fits a fixed number `k` of latent regions to a community matrix or
`CommunityData`. A site can combine several regions: `omega` contains its
mixture weights, while `theta` describes each region's species composition.
Weights express composition, not classification confidence.

```julia
using StableRNGs

gom = fitgom(regional_data; k = 2, nstart = 5, rng = StableRNG(42))
gom.omega                 # sites × regions, each row sums to one
gom.theta                 # species × regions, each column sums to one
gom.membership            # first maximum per site
gom.converged, gom.reason, gom.residual
gom.start_objectives      # all starts, including any that did not converge
gom_indicators = indicators(gom; top_indicators = 3)
```

Counts enter as supplied. Binary matrices fit presence composition; integer
abundances fit a multinomial model. Fractional weights give a weighted
likelihood. The scalar `shape` defaults to `1/(k*species)`; the maximized NEF
objective is `sum(x*log(q)) + shape*sum(log(theta)) + sum(log(omega))/k`, with
`q[i,j] = sum(omega[i,r]*theta[j,r] for r in 1:k)`. `shape` is a positive
pseudocount, not `shape-1`. Rescaling counts changes the prior/data balance.

The optimizer uses EM and accepts log-ratio extrapolations only when they
improve the objective. `accelerate=false` uses EM alone. Each iteration takes
two EM steps. `init_theta` and `init_omega` optionally supply positive factors
for the first start; columns/rows are normalized. Subsequent starts are random.
The largest objective wins. `history`, `start_*` diagnostics, `maxiter`, `tol`
and `accelerate` record the run. `residual` is the largest absolute factor change
under one ordinary EM step; `converged` requires `residual <= tol` (default
`1e-8`). Reasons are `:tolerance`, `:maxiter`, or `:stalled`. Local optima remain
possible. The objective matches scalar-shape maptpx; its conditional SQP/EM
optimizer can follow a different path.

Sites/species keep input order. Regions are sorted by decreasing total
`omega`, with stable exact ties; hard `cluster_sizes` can contain zeroes.
Empty sites are errors. All-zero species columns stay in the model, affecting
the prior and retaining a positive prior contribution in `theta`. CSC input
uses O(nnz + (sites+species)*k) work storage without a dense community or an
nnz × k responsibility array. Inputs are not mutated. Use aggregated
`.community` for lineage or functional fits. `1 <= k <= min(sites,species)`;
automatic k selection, Bayes factors and uncertainty estimates remain in the
long-term plan.

## Regional indicator species

`indicators` ranks species using directed Poisson or Bernoulli KL divergence
between regions. Supply a `GoMResult` or a species-by-region `theta` matrix from an external
model: nonnegative finite intensities for `:poisson`, or probabilities in [0,1]
for `:bernoulli`. Columns are used as supplied, without normalization.

```julia
theta = [0.9 0.1 0.2; 0.2 0.8 0.3; 0.1 0.3 0.7; 0.4 0.4 0.4]
indicator_result = indicators(theta; species_labels = ["A", "B", "C", "D"],
    region_labels = ["west", "east", "south"], method = :bernoulli,
    comparison = :min, top_indicators = 3, shared = false)
indicator_result.top_features  # regions × requested ranks; missing padding
indicator_result.scores
indicator_result.counts        # returned species count per region
```

For target p and comparison q, Poisson KL is `p*log(p/q)+q-p`;
Bernoulli KL adds the corresponding term for `1-p` and `1-q`.
`comparison=:min` (default) uses the least divergence to other regions;
`:max` uses the greatest. At least two regions are required.
With `shared=false`, a species is eligible only for its first maximum-theta
region; `true` lets it appear in several regions. Scores tie in input species
order. Too few eligible species leave `missing` ranks, not repeated species.
`IndicatorSpeciesResult` retains labels, indices, names, scores and options.

Zero/one boundaries use continuous limits, so infinite KL can be valid;
zero against zero gives zero. Scoring uses one species vector per region,
without a species-by-region-by-region tensor. This scorer does not infer theta
from hard partitions; use `fitgom` or supply another fitted statistical model.

## Mean distances across posterior trees

`mean_dist` averages beta distances across posterior trees, aligning sites by
their exact labels. It consumes an iterable once, so a generator can compute
one tree's beta result at a time. All inputs must use the same beta kind and
Sørensen/Jaccard family; `component` selects `:total` (default), `:turnover`, or
`:nestedness`.

```julia
posterior_trees = mktempdir() do dir
    path = joinpath(dir, "posterior.nwk")
    map(["((A:1,B:1):1,(C:1,D:1):1);",
         "((A:1,B:1):2,(C:1,D:1):2);"]) do newick
        write_tree(path, newick)
        read_tree(path)
    end
end

average = mean_dist((phylobeta(regional_data, tree) for tree in posterior_trees))
average.distances             # symmetric mean distances
average.counts                # finite observations per site pair
average.replicates            # number of posterior draws consumed
regions = phyloregion(average; k = 2)
context = prepare_phyloregion(average)
select_linkage(average)
optimal_phyloregion(context)

# Compute and aggregate only these pairs, without full site-by-site matrices.
requested = [(1, 2), (2, 3), (2, 1)]
selected = mean_dist((phylobeta(regional_data, tree; pairs = requested)
    for tree in posterior_trees); pairs = requested)
selected.distances            # vectors preserve request order and repetitions
selected.counts
```

Output order defaults to the first input's labels. Supply `sample_labels` to
set a complete domain explicitly; later labeled inputs may be reordered or
contain subsets, but unknown labels fail. Plain matrices use this common order
and complete domain. For individually labeled distances use named tuples
`(distances=D, sample_labels=labels)`; vectors also need their own `pairs`.
Matrix inputs cannot be mixed with beta results.

By default `nan_policy=:skip` omits NaN separately for each pair. No finite
observations means distance NaN and count zero; this also preserves undefined
diagonals from empty or zero-PD sites. `:propagate` returns NaN if any observed
value is NaN, retaining finite counts. Missing sites/pairs supply no observation.
Infinities, negative distances, nonzero finite diagonals and asymmetric matrices
fail. Different components can have different counts, so their separately
averaged values need not retain the original decomposition identity.

Output `pairs` index the output labels; each input's pair indices refer to its
own labels. Repeated or reversed source pairs must agree and count once per
draw. Pair-vector inputs require explicit output `pairs`. Selected output uses
linear storage in labels and pairs; complete supplied matrices are still
validated. Full output requires quadratic storage. Input arrays are untouched,
and a scaled compensated sum handles very large and very small distances.

Regionalization accepts only complete, finite mean matrices with zero diagonals;
select a defined subset explicitly before clustering. This summarizes distances
across trees and does not estimate posterior support for the resulting regions.
Local R validation covers 42 ensembles, including 18 posterior beta ensembles;
ordinary tests retain four small examples in one reference TSV and do not use R.

## Nonmetric multidimensional scaling

`nmds` converts a complete distance matrix into an ordination using global
NMDS, monotone regression and normalized SMACOF updates. It minimizes Kruskal
stress-1, `sqrt(sum((distance-disparity)^2)/sum(distance^2))`, reported on the
0–1 scale. It also accepts complete beta and `MeanDistResult` inputs.

```julia
using Random

ordination = nmds(beta; component = :total, dims = 2, nstart = 20,
    rng = MersenneTwister(42))
ordination.points             # Sites × dimensions, preserving input row order.
ordination.stress
ordination.converged
ordination.reason
ordination.start_stress       # Final stress of every start.
ordination.stress_history     # Initial and accepted stresses of the best start.

# A PhyloregionResult ordains regions; at least three are needed for dims=2.
regional_ordination = nmds(phyloregion(beta; k = 3);
    nstart = 20, rng = MersenneTwister(42))
```

Beta inputs select `:total`, `:turnover` or `:nestedness`. Functional beta uses
its contained beta result. `nmds(average)` uses the component already selected
by `mean_dist`. `nmds(regions)` uses `region_distances` with labels `"1"`,
`"2"`, ... matching membership IDs, whereas `nmds(context::PhyloregionContext)`
ordains the original sites without rebuilding their clustering.

`nstart` includes the first start: PCoA by default, `init=:random`, or an
explicit sites-by-dimensions coordinate matrix. Subsequent starts are random;
`rng` controls the sequence. The lowest final stress wins. `maxiter=1000`,
`tol=1e-7` (relative stress decrease) and `stress_tol=1e-4` control each run;
`maxiter=0` evaluates starting configurations without updating them. Inspect
per-start iterations, convergence flags and reasons. Numerical convergence
and repeated starts do not establish a global optimum.
Rank-deficient starts can remain in their initial subspace; random starts
explore additional directions.

`ties=:weak` permits equal input dissimilarities to have different fitted
distances, as in vegan's default; `:strong` keeps their fitted values equal.
Only exact ties are grouped. Zero dissimilarities participate in regression
and do not force coincident points. All-equal positive input under weak ties
can give zero stress for arbitrary coordinates. All-zero input is rejected.
Require `1 <= dims < nsites`; an informative NMDS generally needs many more
sites than dimensions. NaN/Inf, asymmetric matrices, nonzero diagonals and
selected-pair results are rejected; resolve undefined sites explicitly.

`NMDSResult` owns centered scores with unit RMS distance to their centroid,
rotated to principal axes. Axis signs, tied principal directions and rotation
are not ecological identifiers. Input dissimilarity units do not set score
units. `dissimilarities`, `distances` and `disparities` contain original,
Euclidean score and monotone fitted distances in upper-triangle column order:
`(1,2), (1,3), (2,3), (1,4), ...`. These aligned vectors support Shepard plots.
No pair-by-iteration history is retained. Pair storage is quadratic; PCoA
initialization adds cubic work, which `init=:random` avoids. Community
transformations, species scores, step-across, half-change scaling and map/color
generation are outside this numeric interface.
