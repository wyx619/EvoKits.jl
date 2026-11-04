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

## 多群落分类 beta

`beta_multi` 利用物种占据数 q 的 `sum(q*(q-1)/2)` 得到全部群落对的共享物种数之和，再通过排序后的群落丰富度累计两两最小值、最大值，扣除共享量得到独有物种数之和。计数和消去采用 `Int128`，最终比值才转为 `Float64`，无需枚举群落对或建立群落间矩阵。

对于 n 个输入群落、p 个物种、m 个所选群落，CSC 路径工作量为 `O(nnz + p + n + m log m)`，辅助存储为 `O(n+p)`；选择子集时仍校验全部输入数值。稠密输入的工作量为 `O(n*p + m log m)`。存储回归在编译后检查 6,000 × 40,000 的 CSC 群落，累计分配上限为 8 MB；这是回归阈值，不是峰值内存或运行时间基准。本轮没有新增计时基准。
