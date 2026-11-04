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
