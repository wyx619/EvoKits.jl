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
`O(s*(e+q)+n^2)` work. Dense preparation scans `n*p` values. Both functions use
`O(q+p+n+e)` auxiliary storage plus the required `O(n^2)` output; no full
tip-distance, tip-pair weight, site-by-branch or dense community matrix is built.
Each target still traverses the full tree. These functions do not currently
accept selected pairs or a prepared context.

Local comparisons with the previous Julia path used balanced unit-branch trees,
eight observed species per site, one Julia/BLAS thread and the median of three
measurements after warm-up:

| Call | Sites × tips | Previous | Current | Previous allocation | Current allocation |
| --- | ---: | ---: | ---: | ---: | ---: |
| Weighted `comdist` | 8 × 1,024 | 354 ms | 0.191 ms | 613 MB | 0.129 MB |
| Weighted `comdistnt` | 24 × 2,048 | 53.3 ms | 0.613 ms | 66.5 MB | 0.166 MB |

These small synthetic cases demonstrate the removed overhead; timings depend
on tree shape, site count and occupancy. Allocations are cumulative, not peak
memory. Separate storage regressions check 12 sites, 20,000 tips and 36 stored
values against a 15 MB allocation bound after compilation for each function.

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
are constructed. For n sites, p tips, e nodes and q observed records, CSC work
is `O(nnz+p+n+e+n*e)`, with `O(q+p+n+e)` auxiliary storage. The full-pool
variance is prepared once in `O(e+p)` work per call. Dense input also scans
`n*p` values. Each site still scans the full tree, and public calls prepare
their inputs again. PSE/PSC and the complete PSD call still use VCV.

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
