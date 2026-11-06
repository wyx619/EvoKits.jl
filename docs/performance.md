# Performance

[Documentation](README.md) · [Home](../README.md) · [简体中文](performance_zh.md)

On this page: [Choose a route](#choose-a-computation-route) · [Measurement](#performance-benchmarks) · [Posterior means](#direct-posterior-beta-means) · [Distances](#community-phylogenetic-distances) · [Moments](#analytical-diversity-moments) · [PD/bootstrap](#pd-variants-and-bootstrap) · [Clustering/NMDS](#regionalization-and-ordination)

Choose a computation route, understand its memory requirements, and interpret
the recorded benchmarks. Models and calling examples live in the topic guides.

## Choose a computation route

| Need | Route | Remaining cost |
| --- | --- | --- |
| Large, sparsely occupied species pool | CSC `CommunityData`; automatic connecting trees | Original input/tree validation and preparation |
| Several PD/PE/beta calls on identical inputs | `prepare_community` | Stored branch sets; output allocation per call |
| Repeated `comdist` / `comdistnt` | `prepare_community_distances`; fix `pairs` if known | One-time preparation and retained workspaces |
| PSV, PSE and PSC together | `psd` | Full-pool variance when enabled |
| Rao summaries without site-pair matrices | `raoD(...; pairwise=false)` | Per-site and pooled calculation |
| Few Expected PD richness values | `expected_pd(...; richness=values)` | Linear tree preparation |
| Posterior-tree uncertainty | `mean_dist(data, trees)` over a `TreeStream` | Per-tree matching/branch sets and requested output |
| A few site comparisons | Supported `pairs` interfaces | Full-input validation; full-community null draws for beta SES |

Dense inputs still require scanning all cells. Complete site-pair outputs have
a quadratic lower bound even when their inputs and workspaces are sparse.

## Runtime

EvoKits loads MKL before other dependencies to select its BLAS backend.
Julia threads and BLAS threads are separate settings. The measurements below
use one of each; tree kernels do not become faster merely by increasing BLAS threads.

## Performance Benchmarks

Tables summarize recorded local measurements, not a fresh benchmark of every
release. Unless stated otherwise: Julia 1.12.7, one Julia/BLAS thread, serial
execution, median of three calls after warm-up. Public-call times include
validation, preparation and output allocation; input construction is excluded.
Some measurements also record an Intel i5-12400. Shape, occupancy, branch scale,
requested output and hardware affect the result; automatic routes need not win
for every input.

Keep these quantities separate:

| Quantity | Meaning |
| --- | --- |
| Cumulative allocation | Bytes allocated during a call, including temporary objects |
| Retained size | Reachable result/context objects measured with `Base.summarysize`; excludes unreferenced inputs |
| Array payload | Analytical byte count, excluding object headers and allocator overhead |
| Process peak RSS | Whole-process maximum, including loading/JIT and earlier work |

MB/GB are decimal. No peak RSS was measured for tables below except the
explicit mean-accumulator process experiment. Historical private-kernel timings
are not public-call benchmarks.

For a small local measurement, start Julia with one thread and run:

```julia
using EvoKits, LinearAlgebra, Statistics
BLAS.set_num_threads(1)
tree = mktempdir() do dir
    path = joinpath(dir, "star.nwk")
    write_tree(path, "(" * join(("sp_$i:1" for i in 1:128), ",") * ");")
    read_tree(path)
end
call() = expected_pd(tree; richness = [10, 50, 100], upper_bound = false)
call()                         # compile and warm up
samples = map(1:3) do _
    GC.gc()                    # outside the timed call
    @timed call()
end
median(s.time for s in samples)
median(s.bytes for s in samples)
Base.summarysize(call())       # retained result, not process peak
```

The notation below uses n sites, p tips, e tree nodes, q positive site–tip
records, r species per site and v working-tree nodes. `nnz` includes stored
zero entries in CSC input. Additional symbols are defined in their sections.

## Posterior tree input

`read_trees` scans once and retains the current record/tree, source paths and
current NEXUS translation table. It does not collect all posterior trees.
`tree_sampler` scans twice, parses only selected records and retains draw
metadata plus sampled trees. Dictionary-based partial Fisher–Yates sampling
avoids an index vector sized to the whole population; repeated selected records
are parsed once and copied.

`CompactTree` conversion builds postorder iteratively in linear work/storage,
preserving child order. Complete `read_tree` calls on actual 65,536-tip left/right
combs took 136/126 ms, allocated 119.73/116.18 MB, and retained 17.89 MB per tree.
This measures reading, parsing and conversion together. Deep-input support does
not imply that every downstream tree operation uses iterative traversal.

## Direct posterior beta means

`mean_dist(data, trees)` prepares community records and requests once, computes
one beta component per tree, and updates a scaled compensated mean. Tree-specific
branch sets and label matching are rebuilt for each draw, with buffer capacity
reused. It avoids a three-component beta result for every tree.

For full outputs with at least 128 sites, sparse branch-to-site memberships
scatter shared lengths into one `n × min(n, block_size)` column buffer. Shared
work follows branch occupancy, `sum(d*(d+1)/2)`, plus quadratic cell work.
Small outputs and selected pairs merge sorted branch sets; near cancellation,
unique branches are summed directly. Deep combs can retain long, highly shared
branch sets, so their complete matrices remain costly.

For b retained site-branch memberships and k requests, selected auxiliary
storage is `O(e+p+q+b+n+k)`. Full output needs
`O(e+p+q+b+n²+n*min(n,block_size))`. A smaller block controls temporary storage,
not the final quadratic mean/count matrices.

### Mean accumulator storage

Both mean entry points retain full sum/count matrices as eventual outputs and
packed lower triangles for scale/correction scratch. Finalization normalizes
and mirrors in place. With 64-bit Int/Float64 the array payload is `24n²+8n`
bytes; returned mean/count arrays alone need `16n²`. At n=6,000 these are
0.864 GB and 0.576 GB respectively. Selected pairs use four one-dimensional
arrays rather than matrix state.

Recorded public calls with 2,048 sites and eight draws:

| Input | Time | Cumulative allocation |
| --- | ---: | ---: |
| Dense distance matrices | 472 ms | 103.03 MB |
| Lazy distance matrices | 586 ms | 103.03 MB |
| Direct posterior beta; 256-tip star, one species/site, block=64 | 236 ms | 102.47 MB |

A separate fresh-process experiment used 4,096 sites, two lazy draws, three
processes and the same load/warm-up sequence. Median cumulative allocation was
404.95 MB and process peak RSS 947.79 MB, including about 551 MB of prior
load/JIT peak. Those totals cannot be inferred from the array-payload formula.

## Community long-table export

For q positive values, `community_longtable` uses `O(nnz+n+p)` work on CSC input
or `O(n*p+n+p)` on dense input. Output and temporary storage are `O(q+n+p)`.
Values are neither densified nor expanded by abundance; full label domains are
copied even for empty inputs. Storage regressions are allocation bounds, not
timing or process-peak measurements.

## Community phylogenetic distances

`comdist` uses postorder/preorder target-weight and distance passes;
`comdistnt` propagates the two nearest distinct target identities. Original
branch segments and stable sibling sums preserve short paths below long stems.
Abundance normalization scales weights first to avoid overflowing sums.

CSC preparation costs `O(nnz+p+n+e)`; dense input additionally scans n*p cells.
Automatic contraction keeps at most `2u-1` nodes for u observed tips, summing
original branches and discarding the common stem above their LCA. Dense
occupancy retains the full tree. The selection policy is an internal heuristic.

After preparation, full output costs `O(t*(v+q)+n²)` for t needed targets.
For k selected requests, sorting/scheduling takes `O(k log k)` and evaluation
`O(t*v+P+k)`, where P sums source richness over unique directional queries.
Full inputs remain validated, but only effective requested endpoints are retained.

### Comdist virtual trees and selected pairs

Recorded full public calls: 100,000 unit-branch tips, 64 CSC sites, three
species/site with abundances `[1,2,3]`, abundance weighting enabled:

| Shape | Time | Cumulative allocation |
| --- | ---: | ---: |
| Balanced | 1.08 ms | 1.74 MB |
| Comb | 1.06 ms | 1.74 MB |
| Star | 0.639 ms | 1.31 MB |

Fewer pairs reduce target passes and output, while direct calls still prepare
the original input. Repeated work benefits from the context below.

### Comdistnt virtual trees and directional requests

The same input with `exclude_conspecifics=true`:

| Shape | Time | Cumulative allocation |
| --- | ---: | ---: |
| Balanced | 1.16 ms | 1.73 MB |
| Comb | 1.21 ms | 1.73 MB |
| Star | 0.662 ms | 1.30 MB |

Each unique pair needs both directions; grouping by target shares tree passes.
Different species at zero distance remain eligible when excluding conspecifics.

### Prepared community distances

`CommunityDistanceContext` retains `O(n+qR+v+k)` records, topology, schedules
and both weighting/distance workspaces; qR is retained presence count. It owns
the snapshot and does not retain the original matrix/tree. Fixed requests cache
schedules and may shrink the topology further. New requests on a complete
context reuse its topology but incur scheduling work. A lock serializes calls
sharing workspaces, and each result owns new arrays.

For the balanced input above, the direct route already compresses the tree:

| Request / function | Direct | Reused | Direct allocation | Reused allocation |
| --- | ---: | ---: | ---: | ---: |
| Full `comdist` | 1.17 ms | 0.235 ms | 1.74 MB | 0.0329 MB |
| Full `comdistnt` | 1.17 ms | 0.227 ms | 1.73 MB | 0.0329 MB |
| Fixed four `comdist` pairs | 0.986 ms | 6.40 μs | 1.64 MB | 96 B |
| Fixed four `comdistnt` pairs | 0.952 ms | 7.70 μs | 1.65 MB | 176 B |

Four requests are `[(1,64),(3,4),(4,3),(1,1)]`. Reused times exclude preparation:
full preparation took 0.978 ms / 1.75 MB allocated and retained 80.0 KB;
fixed-four preparation took 1.13 ms / 1.67 MB and retained 19.8 KB.
Account for these costs when evaluating only one or a few calls.

## Rao diversity partition

RaoD shares sparse records, the distance workspace and automatic connecting
trees. Validation and original-abundance pooling precede compression. With s
nonempty sites, preparation is `O(nnz+p+n+e)` and auxiliary storage
`O(q+p+n+e)`. After preparation, full mode costs `O(s*(v+q)+v+q+n²)`;
summary mode costs `O(s*v+q+v+n)`. Full output needs two n×n matrices;
summary output is linear in n. High occupancy retains full-tree passes.

Recorded abundance-weighted public calls, 32 sites and three species/site,
65,536-tip ultrametric trees, CSC input, summary mode:

| Shape | Time | Cumulative allocation |
| --- | ---: | ---: |
| Balanced | 1.14 ms | 3.215 MB |
| Comb | 1.06 ms | 3.215 MB |
| Star | 0.592 ms | 2.398 MB |

The balanced fixture contracts 131,071 nodes to 191. For 6,000 sites on a
64-tip star, full/summary calls took 579/3.39 ms, allocated 578.02/2.025 MB,
and retained 576.048/0.0482 MB. Tree construction and parsing are excluded.

## Analytical diversity moments

PSV/PSR/PSE/PSC/PSD compute directly from tree summaries, without full VCV,
correlation matrices or per-site submatrices. Centered moments, correlation
deficits and original path segments preserve small contrasts below long stems.
Sparse preparation is linear in `nnz+p+n+e`; dense preparation scans n*p cells.

Expected PD uses subtree sums and probabilities cached by descendant count.
For m requested richness values, preparation is `O(p+e)`, evaluation `O(m*e)`
and storage `O(p+e+m)`. The default full curve is `O(p*e)` work. Empty/all-p
requests validate branches then return in `O(e+m)` work, without tree scratch.
Logarithmic tails preserve representable variance terms after probability
underflow. Exact and upper-bound modes have the same storage order.

Recorded exact variance on 4,096-tip unit-branch trees; compare two available
routes, full curve then select versus `richness=[10,50,100]`:

| Shape | Full curve then select | Requested points | Full allocation | Requested allocation |
| --- | ---: | ---: | ---: | ---: |
| Balanced | 126.22 ms | 0.297 ms | 0.396 MB | 0.297 MB |
| Comb | 344.25 ms | 0.412 ms | 0.686 MB | 0.587 MB |
| Star | 58.02 ms | 0.105 ms | 0.280 MB | 0.182 MB |

The full-curve column was recorded with the preceding full-curve implementation;
default numerical outputs agree exactly. At 65,536 tips the three-point calls
took 4.24/7.91/1.70 ms and allocated 4.72/9.62/2.88 MB for balanced/comb/star
trees. Each three-point result retained 224 B. No large-tree full curve was timed.

### PSV/PSR connecting trees

Sparse site extraction uses a linear-space heavy-path/segment index. In a
v-node working tree, an r-species site costs `O(r log r+r log² v)` to extract
and `O(r)` to evaluate on at most `2r-1` nodes. Dense sites use `O(v+r)`
sequential moments. Full-pool variance is computed once on the **complete
original tree**, including unobserved species, adding `O(e+p)` work/storage.

Recorded public PSV on 100,000 unit-branch tips, 64 CSC sites, three species/site,
correlation scaling enabled:

| Shape | Time without variance | Time with variance | Allocation without/with variance |
| --- | ---: | ---: | ---: |
| Balanced | 5.61 ms | 18.75 ms | 8.19 / 18.59 MB |
| Comb | 6.83 ms | 22.12 ms | 8.19 / 18.59 MB |
| Star | 2.94 ms | 9.29 ms | 4.49 / 10.09 MB |

Variance remains material even with small site trees. Connecting-tree preparation
can increase allocation when variance is enabled. Standalone calls prepare again;
private prepared kernels are not public reuse APIs.

### PSE, PSC and shared PSD

PSE uses weighted cross-subtree contrasts; PSC propagates two distinct target
identities. Contracted paths retain original heights and ordered affine
transforms to preserve correlation and raw-covariance definitions. Rare
underflowing abundance contributions use a BigFloat tree pass.

PSD prepares one site tree for PSV/PSE/PSC, derives PSR from PSV and reuses
PSV for equal-abundance PSE. Extraction/evaluation costs match the preceding
section; optional variance still uses the original species pool.

For the same 100,000-tip/64-site setup, abundances `[1,2,3]`, correlation scaling:

| Shape | PSE time / allocation | PSC time / allocation | PSD with variance time / allocation |
| --- | ---: | ---: | ---: |
| Balanced | 5.69 ms / 8.20 MB | 5.99 ms / 8.20 MB | 20.2 ms / 18.63 MB |
| Comb | 6.62 ms / 8.21 MB | 6.44 ms / 8.20 MB | 22.0 ms / 18.64 MB |
| Star | 2.78 ms / 4.50 MB | 3.51 ms / 4.50 MB | 9.71 ms / 10.13 MB |

Balanced PSD without variance took 6.31 ms / 8.23 MB. Dense collections keep
sequential traversal, so sparse-tree gains should not be extrapolated to full occupancy.

## PD variants and bootstrap

`pd_calc` uses sparse records and descendant counts. Sparse r-tip sites reuse
a linear-space path/site index in `O(r log r+r log² e)` work; dense sites use
`O(e+r)`. TIP directly sums r terminal edges. Preparation is
`O(nnz+p+n+e)` on aligned CSC input; public calls prepare their own index.

Bootstrap scans all p tip probabilities per nonempty draw. A bounded heap keeps
the best r keys or worst p−r keys and their complement, using
`O(min(r,p-r))` heap storage plus selected tips/optional bit markers.
Sampling takes `O(p*(1+log(max(1,min(r,p-r)))))` plus sorting selected keys
on the best-key route. No dense null community is formed; returned replicates
still need `O(n*reps)` storage. RNG consumption and stable tie rules are retained.

Recorded 10,000-tip balanced tree, 16 sites with three species each,
11 bootstrap draws and uniform weights; tree-only calls use one site:

| Call | Time | Cumulative allocation |
| --- | ---: | ---: |
| `pd_calc`, TBL | 0.967 ms | 2.48 MB |
| Community `pd_bootstrap` | 13.0 ms | 2.96 MB |
| Tree `pd_bootstrap` | 2.22 ms | 2.57 MB |

## Null models

Site SES and beta SES retain one null community and reusable species/tree
workspaces, accumulating means and sample SDs online. Specifying beta `pairs`
reduces endpoint branch sets and output, but each draw still randomizes the
complete community to retain the null model's margins and species pool.
Swap effort and runs remain genuine computation costs; storage sparsity does
not make randomization constant-time. See [model definitions](null_models.md).

## Tree distance precision

Weighted distances use root-depth subtraction when relative loss is bounded;
otherwise original branch segments are summed through on-demand weighted
ancestor jumps. Queries take `O(log e)` after indexing; the ancestor/weighted
jump index needs `O(e log e)` storage. This index is separate from the linear
working state used by community kernels. Full tip-pair output remains quadratic
in the number of requested tips.

## Multiple-site taxonomic beta

`beta_multi` derives shared-count sums from species occupancy, then sorts site
richness for the sums of pairwise minima/maxima. Int128 counts postpone floating
conversion until the ratios. No site-pair matrix is created.
For m selected sites, CSC work is `O(nnz+p+n+m log m)`, auxiliary storage
`O(n+p)`; dense work is `O(n*p+m log m)`. Excluded rows remain validated.

## Multiple-site phylogenetic beta

`phylobeta_multi` reuses branch sets and sorted rooted PD. Nonnegative Int128
occupancy/rank coefficients remove shared branches before multiplying lengths.
Close PDs or extreme scales can use a 2304-bit fallback; outputs remain Float64.
For b selected site-branch memberships and m sites, work after branch collection
is `O(b+e+m log m)`, auxiliary storage `O(e+m)`.
`CommunityContext` reuses branch sets. Full input validation still precedes
direct calls; neither route builds dense site-by-branch or site-pair matrices.
The multiple-site storage checks provide allocation bounds, not timing or RSS evidence.

## Regionalization and ordination

The six monotone linkage methods use `O(n²)` nearest-neighbor chaining.
Centroid/median cache row minima with `O(n³)` worst-case time. Clustering
uses a packed Float64 triangle of `n*(n-1)/2` entries; complete inputs and
regional outputs remain quadratic. `PhyloregionContext` retains one triangle
and uses a second temporarily during preparation. Cutting the hierarchy alone
costs `O(n)`; regional summaries still visit site pairs. Linkage evaluation
avoids a full cophenetic matrix. All-cut dispersion scoring visits each pair
once with `O(n)` additional storage.

GoM's CSC working storage is `O(nnz+(n+p)*k)` for k regions, without a dense
community or an nnz × k responsibility array. Indicator scoring uses a species
vector per target region, not a species × regions × regions tensor.

NMDS retains quadratic pair storage, without pair-by-iteration histories.
Default PCoA initialization adds cubic work; `init=:random` avoids that
initialization, while iterations and multiple starts still cost work. Changing
initialization changes the optimization path; inspect stress and convergence.

## Metric value classification

Exact Jenks classification uses multiplicity-weighted dynamic programming over
m distinct sorted values: `O(k*m²)` time and `O(k*m)` workspace for k classes.
`max_work=50_000_000` bounds this workload explicitly. Centering/scaling and
Welford interval variance avoid raw squared-sum cancellation. Equal-width and
quantile methods are alternatives with different statistical definitions;
see [classification rules](conservation.md#metric-value-classes).
