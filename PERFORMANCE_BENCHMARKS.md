# EvoKits Performance Benchmarks

Run date: 2026-07-13

The benchmark uses balanced rooted Newick trees with unit branch lengths and
Julia 1.12.6 on Windows. The first 1k-tip CompactTree timing includes initial
method compilation; later sizes are warm measurements in the same process.

Command shape:

```julia
to_compact_tree(NewickTree.readnw(newick))
_build_tree_query_index(tree)
distances_from_root(tree)
evol_distinct(tree)
sample = zeros(Int, 1, n)
sample[1, 1:div(n, 2)] .= 1
pd(sample, tree)
```

Results:

| Tips | CompactTree | Allocated | Query index | Root distances | ED | Single-site PD |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 0.6037 s* | 51.83 MB | 0.0001 s | 0.0000 s | 0.0000 s | 0.0001 s |
| 10,000 | 0.0169 s | 15.45 MB | 0.0007 s | 0.0000 s | 0.0002 s | 0.0008 s |
| 100,000 | 0.2973 s | 156.42 MB | 0.0133 s | 0.0004 s | 0.0016 s | 0.0136 s |

The 1,000-tip pairwise distance matrix took 0.0086 s and 8.12 MB. It is an
O(n^2) dense matrix and was not allocated for 10,000 or 100,000 tips.

These are smoke benchmarks, not a release performance claim. Repeat them on
the target deployment machine before making allocation or latency promises.
