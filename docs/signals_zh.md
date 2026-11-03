# 系统发育信号

[首页](../README_zh.md) · [English](signals.md)

局部物种池 NRI/NTI、Pagel's lambda、K、独立对比及二元性状 D。

本页示例共用以下四物种树；批量表格示例中的文件路径需替换为实际输入。

```julia
using EvoKits

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
```

## NRI 和 NTI

对于单个观测群落及其本地物种池，使用 `nri_nti`：

```julia
pool = ["A", "B", "C", "D"]
observed = ["A", "B", "D"]

result = nri_nti(tree, pool, observed; nreps = 999, seed = 20260713)
result.nri
result.nti
result.status   # :ok 或 :null_sd_zero
```

零模型从物种池中无放回抽取与观测群落物种数相同的物种。支持 tip 标签和 tip 位置（`1:tree.ntips`）。`nri_nti` 使用稳定随机数生成器和局部虚拟树计算 MPD/MNTD，不假设全球物种池，也不在内存中保存全部零模型重复结果。兼容属性还包括 `NRI`、`NTI`、`MPD_obs` 和 `MNTD_obs`。

批量分析空间网格时，可传入 `grid_id`，使每个网格的随机种子由全局种子稳定派生，不受处理顺序影响。这也是经过验证的 Swq 工作流采用的约定：

```julia
result = nri_nti(tree, pool, observed;
    nreps = 999,
    seed = 20260713,
    grid_id = "1745",
)
```

派生后的实际种子保存在 `result.seed`。不要同时传入 `grid_id` 和显式的 `rng`，因为传入的 RNG 已定义随机数流。

对于大量网格，`nri_nti_table` 接受两张标准分隔长表：一张为本地物种池，一张为观测群落；两表均使用 `grid_id`、`species` 两列。表格可以是 CSV、TSV、TXT 或 CSV.jl 支持的其他分隔格式；需要显式指定制表符时使用 `delim = '\t'`。每个文件中同一网格的记录必须连续。Julia 启动了多个线程时，可并行处理网格。

```julia
batch = nri_nti_table(tree, "local_pool.tsv", "observed_community.tsv";
                        nreps = 999, seed = 20260713, threaded = true,
                        delim = '\t')
batch.grid_ids
batch.results
```

## 连续性状

### Pagel's lambda

`pagel_lambda` 在系统发育 GLS 模型下估计连续单变量性状的 Pagel's lambda。树必须有根、严格二叉且为超度量树；性状值必须有限并按树 tip 顺序提供。

```julia
trait = [1.2, 2.1, 1.8, 2.4]
lambda_result = pagel_lambda(tree, trait)
lambda_result.lambda
lambda_result.ci_lower, lambda_result.ci_upper
lambda_result.p_value
summary_pagel_lambda(lambda_result)
```

估计值限制在 `[0, 1]`：接近零表示考虑拟合残差方差后系统发育协方差较弱，接近一表示保留布朗运动协方差结构。返回的 `PagelLambdaResult` 包含 `success`、`lambda`、`loglik`、`aic`、`sigma2`、`n`、收敛计数、profile likelihood 区间 `ci_lower`/`ci_upper`、似然比检验字段 `lr_stat`/`p_value` 和底层 `fit` 元组。`summary_pagel_lambda` 返回相同信息，并采用 `logLik` 和 `AIC` 字段名。

表格数据可指定性状列，以及可选的分类单元列；分类单元列会在拟合前按 `tree.tip_labels` 对齐：

```julia
using DataFrames

data = DataFrame(species = tree.tip_labels, body_size = trait)
lambda_result = pagel_lambda(
    tree,
    data,
    :body_size;
    taxon_col = :species,
    profile_points = 101,
)
```

`root_model` 接受 `:random` 或 `:fixed`；优化选项包括 `max_iterations`、`rel_tol` 和 `profile_points`。当需要在拟合之外对齐向量或表格时，可使用 `align_traits_to_tree`。

### 系统发育 K（PhyloK）

`phylogenetic_signal_k` 实现 Picante/phytools 的连续性状 K 信号统计量，通过打乱树 tip 上的性状值生成零分布，并与观测 PIC 方差比较：

```julia
using Random

k_result = phylogenetic_signal_k(
    tree,
    trait;
    permutations = 999,
    rng = MersenneTwister(20260717),
)

k_result.k
k_result.p_value
k_result.z_score
k_result.pic_variance_observed
```

返回的 `PhyloSignalResult` 包含 `k`、观测及零模型 PIC 方差汇总、`p_value`、`z_score`、重排后的 `trait`、`tip_labels`、`valid_replicates`、`permutations` 和 `status`。如果性状向量尚未按树排序，可以提供 `tip_labels`：

```julia
k_result = phylogenetic_signal_k(
    tree,
    [1.8, 1.2, 2.1, 2.4];
    tip_labels = ["C", "A", "B", "D"],
    permutations = 999,
)
```

通常，`K` 约为 1 表示接近布朗运动参考模型，`K > 1` 表示更强的系统发育保守性，`K < 1` 表示信号弱于布朗运动预期。p 值根据有效 tip 置换结果的秩计算；显式传入 `AbstractRNG` 可复现置换结果。默认 RNG 为 `StableRNG(20260713)`。

PhyloK 要求有限性状值，以及与 PIC 方差计算相同的严格二叉、正枝长树条件，不接受缺失值。设置 `permutations = 0` 可只计算观测统计量。

独立对比和 Picante/phytools K 信号检验的示例：

```julia
pic = phylogenetic_contrasts(tree, trait)
signal = phylogenetic_signal_k(tree, trait; permutations = 999)
signal.k
signal.p_value
```

PIC 目前要求严格二叉树和正枝长。如需处理多分叉树，请显式使用 `multi2di`；性状缺失值会被拒绝。

## 系统发育 D

`phylo_d` 使用 Fritz-Purvis D 统计量衡量二元性状的系统发育信号，将观测节点变化总和与 tip 状态随机置换、布朗阈值模拟进行比较。

```julia
using Random

trait = [1, 1, 0, 0]  # 按 tree.tip_labels 排序
d = phylo_d(tree, trait; permutations = 999, rng = MersenneTwister(42))

d.d             # D 估计值；0 接近布朗运动，1 接近随机分布
d.p_random      # 随机零模型概率
d.p_brownian    # 布朗运动零模型概率
d.observed
d.mean_random
d.mean_brownian
```

性状必须恰好包含两个有限状态，且不能有缺失值；支持字符串状态和 `Bool`。对于 DataFrame，可使用 `phylo_d(tree, data, :trait; taxon_col = :species)` 按 tip 标签对齐数据。`random_weights` 可为随机零模型指定加权抽样概率。设置 `store_nodal_values = true` 可保留各次置换的节点值矩阵；默认不保留，以降低大型分析的内存占用。兼容属性 `D`、`DEstimate`、`Pval1`、`Pval0`、`Parameters`、`Permutations` 和 `NodalVals` 对应原始 caper 结果名称。

对于支系分析，`phylo_d_subset` 将相同计算应用于经过筛选的内部支系，并为每个支系返回一行汇总。
