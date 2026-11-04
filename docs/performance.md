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
