# 性能

[首页](../README_zh.md) · [English](performance.md)

运行后端与性能基准。

## 运行环境

EvoKits 将 MKL 作为第一个依赖加载，支持选择 Intel MKL 作为 BLAS 后端。

## 性能基准

将 `path` 设置为待测 Newick 文件路径：

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

| Tip 数 | read_tree | 根距离 | ED | 单群落 PD |
| ---: | ---: | ---: | ---: | ---: |
| 10,000 | 0.0170 s | 0.000045 s | 0.000209 s | 0.000807 s |
| 100,000 | 0.2786 s | 0.000302 s | 0.001825 s | 0.013961 s |
| 200,000 | 0.5099 s | 0.059488 s | 0.003887 s | 0.017693 s |
| 500,000 | 1.6090 s | 0.260128 s | 0.010796 s | 0.154415 s |

测量使用单位枝长的平衡有根树，在 Julia 1.12.7 中预热后进行。
