# Performance

[Home](../README.md) · [简体中文](performance_zh.md)

Runtime backend and benchmark timings.

## Runtime

EvoKits loads MKL as its first dependency and supports selecting Intel MKL as
the BLAS backend.

## Performance Benchmarks

Set `path` to the Newick file to measure:

```julia
using EvoKits

tree = read_tree(path)
n = tree.ntips
distances_from_root(tree)
evol_distinct(tree)
sample = zeros(Int, 1, n)
sample[1, 1:div(n, 2)] .= 1
pd(sample, tree)
```

| Tips | read_tree | Root distances | ED | Single-site PD |
| ---: | ---: | ---: | ---: | ---: |
| 10,000 | 0.0170 s | 0.000045 s | 0.000209 s | 0.000807 s |
| 100,000 | 0.2786 s | 0.000302 s | 0.001825 s | 0.013961 s |
| 200,000 | 0.5099 s | 0.059488 s | 0.003887 s | 0.017693 s |
| 500,000 | 1.6090 s | 0.260128 s | 0.010796 s | 0.154415 s |

Measurements use balanced rooted trees with unit branch lengths on Julia 1.12.7
after a warm-up.

## Posterior tree input

`read_trees` scans once, retaining the current command/record and parsed tree
plus source paths and the current NEXUS translation table. Metadata storage
does not grow with the number of trees; skipped blocks are framed command by
command. `tree_sampler` scans input twice, parses only selected records,
and retains O(n) draw metadata plus sampled trees. Without-replacement indices
use a partial Fisher-Yates shuffle over a dictionary, with no array sized to
the posterior population. Repeated draws are parsed once and copied.

Storage regressions check a stream after 2,000 records, three sampled trees
from 5,000 records in Newick and NEXUS, and eight distinct indices from a population of 10^12.
These bound retained objects/index allocations, not peak process memory or
runtime. The downstream `mean_dist` generator retains no collection of beta
matrices; full site distances still require quadratic storage, while selected
pairs use linear storage in sites and requested pairs.

## Community long-table export

For n sites, p species and q positive records, `community_longtable` takes
`O(nnz+n+p)` work for CSC input and `O(n*p+n+p)` for dense input, including
validation. Output and temporary storage are `O(q+n+p)`. CSC input is scanned
directly; neither a dense matrix nor an abundance-expanded occurrence table
is created. Label domains are copied even if there are no positive records.

The storage regression exports four values from a 10,000 × 50,000 CSC matrix,
with one abundance of 1e100, against a 40 MB cumulative-allocation bound after
compilation and a 6 MB retained-output bound. These are regression limits,
not peak-memory or timing measurements.

## Community phylogenetic distances

`comdist` processes one target community at a time. A postorder pass sums
target weights and weighted distances below each node; a preorder pass adds
contributions from outside each subtree. Prefix/suffix sibling sums avoid
subtraction of nearly equal masses and long shared root paths. `comdistnt`
propagates the two nearest distinct target identities in two passes, retaining
an eligible alternative when excluding the source species. Row weights are
scaled before normalization to avoid overflowing abundance sums/products.

For n sites, s nonempty targets, p tips, e tree nodes and q observed site–tip
records, CSC preparation takes `O(nnz+p+n+e)` work, followed by
`O(s*(e+q)+n^2)` work for the full-tree route. Dense preparation scans `n*p`
values. No full tip-distance, tip-pair weight, site-by-branch or dense community
matrix is built. Both functions can reuse a compressed connecting tree within
the call. `prepare_community_distances` also retains prepared records,
topology and workspaces for repeated calls.

Both functions can contract the tree once to their u observed tips and connecting
nodes, retaining at most `2u-1` nodes. A linear postorder pass removes
unobserved subtrees and contracts unary paths by adding original edge lengths.
The common stem above the observed LCA is discarded. This gives cancellation
safe segments on deep combs without constructing a binary LCA index or sorting
tip IDs. Each function's two-pass kernel then visits only the compressed tree.

With `pairs`, validation still covers the full input, while row preparation
retains requested endpoints only. Unordered pairs are deduplicated;
`comdist` groups canonical targets, while `comdistnt` schedules both directional
queries by target. Each needed target is traversed once and each unique pair
is evaluated once; reversed/repeated values are identical. Endpoints
paired only with empty sites do not enlarge the connecting tree. Let k be
the number of requests, t the number of distinct valid canonical targets,
qR the retained records and v the working tree size. After preparation and
`O(e+qR)` contraction, selected work is `O(k log k+t*v+P+k)`, where P sums
the source-site richness over unique directional queries (`comdistnt` visits
both directions). Auxiliary storage is
`O(p+n+e+qR+k+v)`, and output is a length-k vector. Full output still needs
`O(n^2)` storage and `O(t*(v+q)+n^2)` work after preparation.

The automatic route uses compression for trees of at least 256 nodes when
`2u-1 <= e/4`, tightened to `e/16` for a single target. Serial calibration
compares both routes over tree shapes, site counts and occupancy; dense
occupancy stays on the full-tree route to avoid compression allocations.
Occupancy checking stops once the compression threshold has been exceeded,
avoiding scans of remaining repeated records on dense inputs.
The policy is conservative and does not claim the fastest route for every
dataset. Each direct matrix/data call still scans the original tree/input;
prepared-context calls reuse that work.

Earlier dense-to-sparse comparisons used balanced unit-branch trees,
eight observed species per site, one Julia/BLAS thread and the median of three
measurements after warm-up:

| Call | Sites × tips | Dense path | Sparse tree passes | Dense allocation | Sparse allocation |
| --- | ---: | ---: | ---: | ---: | ---: |
| Weighted `comdist` | 8 × 1,024 | 354 ms | 0.191 ms | 613 MB | 0.129 MB |
| Weighted `comdistnt` | 24 × 2,048 | 53.3 ms | 0.613 ms | 66.5 MB | 0.166 MB |

These small synthetic cases demonstrate the removed overhead; timings depend
on tree shape, site count and occupancy. Allocations are cumulative, not peak
memory. Separate storage regressions check 12 sites, 20,000 tips and 36 stored
values against a 15 MB allocation bound after compilation for each function.
Selected-pair storage also checks 6,000 sites × 40,000 tips with 18,000 stored
values: four requests retain 12 endpoint records and a 13-node star, under a
10 MB cumulative-allocation bound. These are regression bounds, not peak memory.

### Comdist virtual trees and selected pairs

The following measurements compare the preceding sparse full-tree implementation
with automatic compression. They use Julia 1.12.7 on an i5-12400, one Julia/BLAS
thread and the median of three calls after warm-up. Trees have unit branches;
each site has three observations with abundances `[1,2,3]`, and abundance
weighting is enabled. Tree construction is excluded; public-call preparation
and the complete matrix output are included:

| Shape | Tips | Sites | Full-tree baseline | Automatic | Baseline allocation | Automatic allocation |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Balanced | 100,000 | 8 | 36.1 ms | 1.15 ms | 12.0 MB | 1.64 MB |
| Balanced | 100,000 | 64 | 300 ms | 1.08 ms | 12.1 MB | 1.74 MB |
| Comb | 100,000 | 8 | 40.1 ms | 0.782 ms | 12.0 MB | 1.64 MB |
| Comb | 100,000 | 64 | 317 ms | 1.06 ms | 12.1 MB | 1.74 MB |
| Star | 100,000 | 8 | 14.5 ms | 0.518 ms | 6.40 MB | 1.22 MB |
| Star | 100,000 | 64 | 106 ms | 0.639 ms | 6.45 MB | 1.31 MB |

Private route comparisons on the 100,000-tip balanced tree with 64 sites
give the following selected-pair measurements. Input preparation is included;
pair-argument parsing is outside timing. Four requests are
`[(1,64),(3,4),(4,3),(1,1)]`; the chain requests adjacent site indices:

| Requests | Forced full tree | Automatic | Full-tree allocation | Automatic allocation |
| --- | ---: | ---: | ---: | ---: |
| One pair | 7.60 ms | 1.06 ms | 12.0 MB | 1.64 MB |
| Four requests | 13.9 ms | 1.03 ms | 12.0 MB | 1.64 MB |
| Chain, 63 pairs | 306 ms | 1.14 ms | 12.0 MB | 1.71 MB |

Measured separately on the same balanced input, record preparation takes
0.330 ms / 0.820 MB and contraction to 383 nodes takes 0.603 ms / 0.842 MB.
With already allocated workspaces, one target pass takes 4.67 ms on the full
tree and 0.00640 ms on the compressed tree; both allocate zero bytes. These
internal passes exclude preparation; the public context described below
retains the complete preparation and workspaces.

Direct matrix/data calls still require linear preparation, and complete
matrices retain quadratic output storage. Denser cases
retain the full route: a 1,000-tip star with 64 sites changed from 0.521 to
0.535 ms. Timings vary with occupancy, shape and hardware; the policy does
not guarantee an improvement for every input. All MB figures are decimal
cumulative allocations, not peak memory.

### Comdistnt virtual trees and directional requests

The same serial protocol compares `comdistnt` with its preceding sparse
full-tree implementation. Both exclusion settings are measured; this table
uses abundance weighting and `exclude_conspecifics=true`, three observations
per site and complete public-call matrix output, including preparation:

| Shape | Tips | Sites | Full-tree baseline | Automatic | Baseline allocation | Automatic allocation |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Balanced | 100,000 | 8 | 9.98 ms | 1.00 ms | 7.20 MB | 1.64 MB |
| Balanced | 100,000 | 64 | 68.9 ms | 1.16 ms | 7.25 MB | 1.73 MB |
| Comb | 100,000 | 8 | 10.1 ms | 0.830 ms | 7.20 MB | 1.64 MB |
| Comb | 100,000 | 64 | 72.7 ms | 1.21 ms | 7.25 MB | 1.73 MB |
| Star | 100,000 | 8 | 4.67 ms | 0.489 ms | 4.00 MB | 1.22 MB |
| Star | 100,000 | 64 | 34.4 ms | 0.662 ms | 4.05 MB | 1.30 MB |

For the same 100,000-tip balanced input with 64 sites, selected requests
include input preparation but exclude pair-argument parsing. The four-request
and chain definitions match the preceding `comdist` table. Each unique pair
requires both directions; target passes are shared between requests:

| Requests | Forced full tree | Automatic | Full-tree allocation | Automatic allocation |
| --- | ---: | ---: | ---: | ---: |
| One pair | 3.94 ms | 0.978 ms | 7.21 MB | 1.64 MB |
| Four requests | 5.75 ms | 1.08 ms | 7.21 MB | 1.65 MB |
| Chain, 63 pairs | 69.0 ms | 1.09 ms | 7.23 MB | 1.71 MB |

Measured independently on the complete balanced input, record preparation
takes 0.329 ms / 0.820 MB and contraction to 383 nodes takes
0.600 ms / 0.842 MB. With prepared workspaces, one target pass takes
1.21 ms on the full tree and 0.00480 ms on the compressed tree, both with
zero allocated bytes. These passes exclude preparation and directional
source averaging; the public context reuses preparation and workspaces.

Sparse large-tree improvements do not imply a speedup for every input.
The 1,000-tip star with 64 sites stays on the full route; without exclusion
it changed from 0.448 to 0.470 ms. Full output still needs quadratic storage,
and each direct matrix/data call still validates/prepares the complete original input.
All allocations are cumulative, not peak memory.

### Prepared community distances

`prepare_community_distances` snapshots aligned abundance/presence records,
one working topology and both distance workspaces. It retains
`O(n+qR+v+k)` storage for sites, retained observations, working-tree nodes
and fixed requests. The original matrix/tree and their full indices are not
retained. Dense collections may keep the full topology; sparse collections
can keep a compressed one. Unlike `CommunityContext`, no per-site branch
sets or occupancy summaries are prepared for these distances.

Complete repeated calls cost `O(t*(v+qR)+n^2)` including new output. Fixed-pair
calls cost `O(t*v+P+k)`, with cached ordering and directional schedules.
Changing pair requests on a complete context still costs `O(k log k+n)` to
schedule, but reuses its records, topology and workspaces. Such requests use
the topology for the complete prepared collection; preparing fixed pairs
can reduce the working tree further when only a few endpoints are needed.

Preparation validates the full input once and owns copied values/topology.
Repeated calls validate only new pair arguments. A lock serializes calls
sharing the workspaces, and every result array is independent. Complete
outputs still require quadratic storage, and preparing both weighting modes
and workspaces costs time and retained memory up front.

The following public-call comparisons use the same serial protocol,
100,000 tips, 64 sites and three observations per site with abundances
`[1,2,3]`. Weighting is enabled; nearest distance also excludes conspecifics.
The direct baseline already includes automatic tree compression. Reused
times exclude the one-time preparation:

| Shape | Call | Direct | Reused | Direct allocation | Reused allocation |
| --- | --- | ---: | ---: | ---: | ---: |
| Balanced | `comdist` | 1.17 ms | 0.235 ms | 1.74 MB | 0.0329 MB |
| Balanced | `comdistnt` | 1.17 ms | 0.227 ms | 1.73 MB | 0.0329 MB |
| Comb | `comdist` | 1.19 ms | 0.256 ms | 1.74 MB | 0.0329 MB |
| Comb | `comdistnt` | 1.10 ms | 0.260 ms | 1.73 MB | 0.0329 MB |
| Star | `comdist` | 0.622 ms | 0.0966 ms | 1.31 MB | 0.0329 MB |
| Star | `comdistnt` | 0.667 ms | 0.156 ms | 1.30 MB | 0.0329 MB |

On the balanced input, fixing four requests during preparation retains
23 working-tree nodes. The requests are `[(1,64),(3,4),(4,3),(1,1)]`;
the chain connects all adjacent site indices:

| Request | Call | Direct | Reused | Direct allocation | Reused allocation |
| --- | --- | ---: | ---: | ---: | ---: |
| Four | `comdist` | 0.986 ms | 6.40 μs | 1.64 MB | 96 B |
| Four | `comdistnt` | 0.952 ms | 7.70 μs | 1.65 MB | 176 B |
| Chain, 63 | `comdist` | 1.16 ms | 192 μs | 1.71 MB | 576 B |
| Chain, 63 | `comdistnt` | 1.11 ms | 139 μs | 1.71 MB | 1,152 B |

Changing a single pair on the complete prepared context, whose topology
has 383 nodes, takes 10.7 μs / 240 B for `comdist` and
12.8 μs / 640 B for `comdistnt`. This still reuses preparation but builds
fresh request schedules. A separately prepared fixed single pair retains
11 nodes instead.

Measured separately on that balanced input, complete preparation takes
0.978 ms / 1.75 MB cumulative allocation and retains 80.0 KB, including
both distance workspaces and weighting modes. Fixed-four preparation takes
1.13 ms / 1.67 MB and retains 19.8 KB. These are one-time preparation costs;
retained object sizes exclude the original input/tree and are not peak memory.

Dense cases retain most tree work. On a fully occupied 4,096-tip balanced
tree with eight sites, `comdist` changes from 1.48 to 0.857 ms and
`comdistnt` from 2.34 to 1.80 ms; reused allocation is 592 B each. Preparing
that complete dense context takes 2.87 ms / 4.93 MB and retains 2.59 MB.
All figures measure cumulative allocations, not peak memory. Reuse is
intended for repeated analysis; preparation and retained storage must also
be counted when evaluating one or a few calls.

## Rao diversity partition

`raoD` reuses the sparse row preparation and tree workspace of `comdist` for
`Dkl`, takes its diagonal for `Dkk`, and evaluates the pooled community once
with the same two-pass kernel. Original abundances are scaled before pooling, so
finite values can be used even when their raw sum would overflow. No dense
community, full tip-distance or tip-pair weight matrix is constructed.

For n sites, s nonempty sites, p tips, e tree nodes and q observed records,
CSC preparation takes `O(nnz+p+n+e)` work, followed by
`O(s*(e+q)+e+q+n^2)`. Auxiliary storage is `O(q+p+n+e)`; the public `Dkl`
and `H` still require two `O(n^2)` matrices. Dense preparation scans `n*p`
values. Each target and the pooled community traverse the full tree.

Local comparisons with the previous Julia RaoD implementation used balanced
unit-branch ultrametric trees, eight observed species per site, one Julia/BLAS
thread and the median of three measurements after warm-up:

| Weighted `raoD` | Previous | Current | Previous allocation | Current allocation |
| --- | ---: | ---: | ---: | ---: |
| 12 sites × 256 tips | 50.6 ms | 0.107 ms | 166 MB | 0.0512 MB |
| 8 sites × 1,024 tips | 958 ms | 0.231 ms | 1,235 MB | 0.164 MB |

These synthetic examples measure cumulative allocations, not peak memory;
results depend on tree shape, site count and occupancy. A separate storage
regression checks 12 sites, 20,000 tips and 36 stored values against a 15 MB
allocation bound after compilation in both weighting modes.

## Analytical diversity moments

`psv` and `psr`, including `compute_var=true`, use subtree pair moments and
rerooted tip means. Centered means/SDs avoid subtracting large nearly equal
moments; scaled covariance uses correlation deficits to retain tiny contrasts
below long stems. Neither full VCV/correlation matrices nor site submatrices
are constructed by PSV/PSR/PSE/PSC/PSD. For n sites, p tips, e nodes and q
observed records, aligned CSC work
is `O(nnz+p+n+e+n*e)`, with `O(q+p+n+e)` auxiliary storage. The full-pool
variance is prepared once in `O(e+p)` work per call. Dense input also scans
`n*p` values. Each site still scans the full tree, and public calls prepare
their inputs again. Permuted labels additionally sort per-site presence
records for PSV/PSR/PSC; weighted PSE/PSD preparation needs no such sort.

`pse` accumulates weighted cross-subtree contrasts instead of subtracting a
quadratic form from a weighted trace. Abundances are scaled before summation;
if positive normalized weights underflow, a rare BigFloat tree pass preserves
their contributions. `psc` propagates the two nearest distinct identities
under `1-Cij` in two tree passes. Its raw-covariance path uses shared node
heights, preserving the original definition even when the value is negative.

`psd` shares one label/record preparation, tree geometry, full-pool variance
and reusable mean/mass workspace. It derives PSR and its variance from PSV,
and reuses PSV for equal-abundance PSE. Output contains seven independent
vectors in site order; it needs no site-pair matrix. Separate public calls
still prepare their inputs independently.

`expected_pd` / `pd_variance` replace ancestor–descendant edge-pair enumeration
with subtree sums. The full Bernoulli richness curve takes `O(p*e)` work and
`O(p+e)` storage, replacing `O(p*e^2)` work. Exact variance and the covariance
upper bound use the same storage order. Probabilities are cached by descendant
count; logarithmic tail evaluation retains representable variance terms when
absence probabilities underflow. The distinct sampling models and corrected
PSV variance are described in the [community guide](community.md#psv-psr-and-analytical-pd-moments).

Local comparisons use balanced unit-branch trees, Julia 1.12.7, one Julia/BLAS
thread and the median of three calls after warm-up. PSV uses eight sites with
three presences each, scaled covariance and variance enabled; expected PD
uses the default upper bound. Tree construction is excluded:

| Call | Tips | Previous | Current | Previous allocation | Current allocation |
| --- | ---: | ---: | ---: | ---: | ---: |
| `expected_pd` | 256 | 84.2 ms | 0.451 ms | 1.11 MB | 0.0258 MB |
| `psv`, including variance | 512 | 399 ms | 0.151 ms | 559 MB | 0.0924 MB |

The following current-path measurements use exact expected PD variance and
the same PSV settings. Balanced, comb and star topologies are constructed
iteratively outside timing; the benchmark also varies site richness from
3 to 300 and disables PSV variance in separate cases.

| Call | Shape | Tips | Sites | Time | Allocation |
| --- | --- | ---: | ---: | ---: | ---: |
| `expected_pd`, exact | Balanced | 256 | — | 0.560 ms | 0.0259 MB |
| `expected_pd`, exact | Balanced | 1,024 | — | 10.5 ms | 0.100 MB |
| `expected_pd`, exact | Balanced | 4,096 | — | 140 ms | 0.396 MB |
| `expected_pd`, exact | Comb | 4,096 | — | 390 ms | 0.686 MB |
| `expected_pd`, exact | Star | 4,096 | — | 64.4 ms | 0.280 MB |
| `psv`, including variance | Balanced | 1,000 | 8 | 0.197 ms | 0.178 MB |
| `psv`, including variance | Balanced | 10,000 | 8 | 1.95 ms | 1.76 MB |
| `psv`, including variance | Balanced | 100,000 | 8 | 39.7 ms | 17.6 MB |
| `psv`, including variance | Comb | 100,000 | 8 | 37.0 ms | 17.6 MB |
| `psv`, including variance | Star | 100,000 | 8 | 18.2 ms | 9.60 MB |

For the 100,000-tip balanced tree, preparing scaled tree geometry alone takes
2.80 ms / 6.40 MB. Repeating the internal single-site mean kernel with that
geometry and allocated workspace takes 2.88 ms / zero allocated bytes. This
kernel excludes label preparation and full-pool variance; it is not a public
prepared-context interface. All allocations are cumulative, not peak memory;
timings depend on tree shape, sites, occupancy and hardware. A full 100,000-tip
Float64 VCV alone would require 80 GB. Separate numerical/storage regressions
use exhaustive small-tree subsets, extreme branch scales and a 20,000-tip
sparse PSV input.

### PSE, PSC and shared PSD

The following comparisons use the same serial protocol, balanced unit-branch
trees, eight sites and three presences per site with abundances `[1,2,3]`.
Correlation scaling is enabled, and PSD includes variance:

| Call | Tips | Previous | Current | Previous allocation | Current allocation |
| --- | ---: | ---: | ---: | ---: | ---: |
| `pse` | 1,024 | 13.3 ms | 0.109 ms | 25.4 MB | 0.111 MB |
| `psc` | 1,024 | 13.3 ms | 0.162 ms | 25.4 MB | 0.125 MB |
| `psd`, including variance | 1,024 | 28.6 ms | 0.384 ms | 51.3 MB | 0.235 MB |

Current calls with the same site settings scale as follows. Separate cases
also cover one site, raw covariance and richness 30/300:

| Shape | Tips | PSE | PSC | PSD, including variance | PSD allocation |
| --- | ---: | ---: | ---: | ---: | ---: |
| Balanced | 1,000 | 0.105 ms | 0.177 ms | 0.371 ms | 0.230 MB |
| Balanced | 10,000 | 0.980 ms | 1.24 ms | 3.39 ms | 2.25 MB |
| Balanced | 100,000 | 21.0 ms | 14.9 ms | 65.9 ms | 22.4 MB |
| Comb | 100,000 | 20.7 ms | 19.4 ms | 59.9 ms | 22.4 MB |
| Star | 100,000 | 8.86 ms | 6.95 ms | 28.6 ms | 12.0 MB |

On the 100,000-tip balanced tree, four separate calls to PSV/PSR/PSE/PSC
take 107.7 ms / 57.6 MB, versus 65.4 ms / 22.4 MB for PSD. With equal
abundances, PSD takes 51.8 ms because PSE reuses PSV. Preparing the weighted
records alone takes 0.215 ms / 0.803 MB; tree geometry takes 2.52 ms / 6.40 MB.
With prepared geometry and workspace, the internal single-site PSE and PSC
kernels take 2.30 ms and 1.59 ms, respectively, with zero allocated bytes.
These internal calls exclude input preparation and full-pool variance and
do not expose a public prepared context. All figures are cumulative
allocations, not peak memory. Storage regressions cover both scaling modes
at 20,000 tips, including PSD's allocation savings over separate calls.

## Tree distance precision

Weighted tree distances use root-depth differences only when relative loss is
bounded. Otherwise they sum original branch segments through weighted ancestor
jumps, built on demand and reused within matrix/vector calls. Path queries
remain `O(log e)` after indexing, including on deep combs. The ancestor index
and optional weighted jumps use `O(e log e)` storage; full pairwise output
still requires quadratic space in the number of requested tips.

## Multiple-site taxonomic beta

`beta_multi` derives pairwise shared counts from species occupancy `q` as
`sum(q*(q-1)/2)`. Sorting selected-site richness supplies the sums of pairwise
minima and maxima; subtracting shared counts gives the required unique-count
sums. Counts and cancellation use `Int128`, with conversion to `Float64` at
the final ratios. Neither pair enumeration nor a site-pair matrix is needed.

For n input sites, p species and m selected sites, CSC work is
`O(nnz + p + n + m log m)` and auxiliary storage is `O(n+p)`; the full input is
validated even for a subset. Dense input takes `O(n*p + m log m)` work.
The storage regression checks a 6,000 × 40,000 CSC collection against an
8 MB cumulative-allocation bound after compilation. This is a regression bound,
not peak-memory or runtime benchmark evidence. No new timing benchmark was run.

## Multiple-site phylogenetic beta

`phylobeta_multi` reuses sorted branch sets and orders selected sites by rooted
PD. For an edge occupied at q sites, let R be the sum of zero-based site ranks
where the edge occurs. Its coefficients in maximum/minimum unique-length sums
are `R-q*(q-1)/2` and `q*(m-1)-R-q*(q-1)/2`. These nonnegative `Int128`
coefficients remove shared branches before multiplying by branch length,
preserving small terminal lengths attached to long stems. Adjacent PD gaps
supply the difference between the two unique-length sums, and nestedness uses
the defining product. Close PDs, overflowing intermediate sums or extreme
scales use a 2304-bit arithmetic fallback; outputs remain `Float64`.

After branch collection, m selected sites, e tree edges and b site-branch
memberships require `O(b + e + m log m)` work and `O(e+m)` auxiliary storage.
Matrix input first validates all rows and collects presence/branches;
`CommunityContext` reuses its branch sets and tree snapshot. Neither path
builds dense community, site-by-branch or site-pair matrices. The storage
regression checks a 6,000 × 4,000 CSC collection with 18,000 stored presences,
against 20 MB cumulative allocations for a direct call and 5 MB for a prepared
call after compilation. These are regression bounds, not measured peak memory
or timing comparisons.
