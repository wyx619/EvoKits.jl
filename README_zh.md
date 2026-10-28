# EvoKits.jl

EvoKits 是一个原生 Julia 系统发育分析工具包，提供树的读写、紧凑的内存树表示、拓扑操作、群落系统发育指标、树节点与边的身份映射，以及轻量级比较分析功能。

## 安装与加载

```julia
using Pkg
Pkg.add(url = "https://github.com/wyx619/EvoKits.jl")
using EvoKits
```

EvoKits 将 MKL 作为第一个依赖加载，支持选择 Intel MKL 作为 BLAS 后端。

## 树的读写与拓扑

```julia
using EvoKits

tree = read_tree("example.nwk")

ultrametric_tree = make_ultrametric(tree)
binary_tree = multi2di(ultrametric_tree; random = false)
write_tree("binary.nwk", binary_tree)
```

`CompactTree` 是 EvoKits 工具的通用输入类型。它缓存拓扑、枝长、遍历顺序和 tip 标签查找信息。`read_tree` 将 Newick 文件直接读取为这种表示。

树查询明确使用紧凑节点 ID、tip 位置或 tip 标签：

```julia
mrca(tree, "A", "B")
distance(tree, 1, 2)
pairwise_distance(tree)
descendant_tips(tree, tree.root; labels = true)
tree_imbalance(tree; metric = :sackin)
```

加权距离使用枝长；设置 `weighted = false` 可改为计算边数。源代码分为 `base`、`tree`、`traits` 和 `community` 四层。

## 树的修剪与映射

```julia
pruned = keep_tip(tree, ["A", "C", "D"])
without_b = drop_tip(tree, "B")

map = build_phyloref(pruned)
node_map = phyloref_node_table(pruned; map = map)
edge_map = phyloref_edge_table(pruned; map = map)

ape_edges = ape_edge_matrix(pruned; map = map)
signature = phylo_edge_signature(pruned, 1; map = map)
```

`PhyloRef` 提供稳定的节点和边身份标识，包括 `CompactTree` ID 与 R/ape 编号之间的转换。映射表和签名方便在不同工具及运行环境之间比较树结果。

## 群落指标

```julia
community = [
    1 0 1
    0 1 1
]

pd_result = pd(community, tree; tip_labels = tree.tip_labels)
ses_result = ses_mpd(community, tree; runs = 999)
```

群落工具包括系统发育多样性 `pd`、进化独特性 `evol_distinct`，以及 PD、MPD 和 MNTD 的标准化效应量指标 `ses_pd`、`ses_mpd` 和 `ses_mntd`。

群落 beta 多样性及其他多样性指标使用相同的行列约定：

```julia
comdist(community, tree)
comdistnt(community, tree)
unifrac(community, tree)
phylosor(community, tree)
raoD(community, tree)

psv_result = psv(community, tree)
psv_result.psv
expected_pd(tree).expected_pd
```

`raoD` 返回 `RaoDResult`，包含 `Dkk`、`Dkl`、`H`、`alpha`、`beta` 和 `Fst`，遵循 Picante 的分解约定。它要求超度量树，并使用到最近共同祖先的距离（`patristic / 2`）。

矩阵的行表示群落，列表示树的 tip。如果输入列的顺序与树 tip 顺序不同，请提供 `tip_labels`。当距离或多样性未定义时，空群落返回 `NaN`；PD、UniFrac、PhyloSor 和 RaoD 接受单物种群落。

非丰度加权的 MPD/MNTD 使用虚拟树计算，无需为每个群落扫描大树的所有边。零模型均值和标准差使用 Welford 在线算法累计。

如需 Caper 兼容的多分叉 ED 处理方式，可使用 `evol_distinct(tree; type = :fair_proportion, polytomy_correction = :isaac)`；其他支持的模式为 `:mooers` 和 `:none`。

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

## 比较分析工具

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

## EDGE2

EDGE2 将进化独特性与物种灭绝概率结合，使用与其他工具相同的 `CompactTree`，无需额外的树解析器或 `Phylo.jl` 依赖。

```julia
pext = [0.10, 0.20, 0.30, 0.40]  # 按 tree.tip_labels 排序
edge2 = edge2_mod(tree, pext)

edge2.species
edge2.edge       # EDGE2 得分
edge2.ed         # 进化独特性
edge2.pd         # 系统发育多样性总量
edge2.epd_loss   # 预期 PD 损失
```

参考实现中的字段拼写 `EDGE`、`ED`、`TBL`、`PD` 和 `ePDloss` 仍可作为兼容属性使用。灭绝概率必须是 `(0, 1]` 内的有限值，并按树 tip 顺序排列。

基于类别的工作流可使用 `PextDistribution`、`calc_pext_distribution`、`load_pext_distribution_csv` 和 `sample_pext`。`reorder_pext_to_tree` 会先验证物种覆盖情况，再按树顺序返回概率。

## 完整示例

以下示例在同一会话中使用树、群落和性状分析功能：

```julia
using EvoKits
using Random

write_tree("example.nwk", "((A:1,B:1):1,(C:1,D:1):1);")
tree = read_tree("example.nwk")

community = [
    1 1 0 0
    0 1 1 0
]
pd_result = pd(community, tree)
beta = phylosor(community, tree)
rao = raoD(community, tree)

trait = [1.0, 2.0, 3.0, 4.0]
pic = phylogenetic_contrasts(tree, trait)
signal = phylogenetic_signal_k(tree, trait;
    permutations = 99, rng = MersenneTwister(7))

pd_result.pd
beta[1, 2]
rao.Dkk
pic.contrasts
signal.k
```

## 公共 API 约定

### 共同排序约定

`tree.tip_labels` 是标准 tip 顺序。tip 位置为 `1:tree.ntips` 内的整数；紧凑节点 ID 为 `1:tree.nnodes` 内的整数。接受标签的函数通过 `tree.tipname_to_id` 查找。

群落矩阵以行为群落、列为 tip。所有返回向量或矩阵均保持群落行顺序。`tip_labels` 描述输入列顺序，并将列重排至 `tree.tip_labels` 顺序；它不会重命名或重排群落行。重复、缺失或未知的标签会触发 `ArgumentError`。

群落值必须有限且非负。大于零表示存在；丰度加权方法将正值作为权重。

性状向量默认按树顺序排列。提供 `tip_labels` 时，它描述所提供向量的顺序，程序会将向量重排至 `tree.tip_labels` 顺序。K 和 PIC 拒绝缺失值及非有限值。

### 矩阵语义

`comdist`、`comdistnt`、`unifrac` 和 `phylosor` 返回普通的 `n x n` `Float64` 矩阵，索引顺序与群落行顺序一致。

- `comdist`：两个群落之间的平均成对距离；对角线为群落内平均距离；空群落返回 `NaN`。
- `comdistnt`：两个方向最近分类单元平均距离的对称平均值；对角线为群落内最近分类单元平均距离。设置 `exclude_conspecifics=true` 时，单物种群落的对角线为 `NaN`。
- `unifrac`：非加权枝长距离；对角线为零，包括空群落之间的比较。当非空群落跨越正枝长时，空群落与非空群落之间的距离为一。
- `phylosor`：`2 * shared_PD / (PD1 + PD2)`；非空群落的对角线为一，空群落为零。

`raoD` 返回包含 `Dkk`、`Dkl`、`H`、`total`、`alpha`、`beta` 和 `Fst` 的 `RaoDResult`，遵循 Picante 分解约定。它要求超度量树，使用到最近共同祖先的距离（`patristic / 2`）。

### 群落结果类型

- `pd` 返回 `PDResult`，包含 `sample_labels`、`pd`、`richness` 和 `include_root`。
- `ses_pd`、`ses_mpd` 和 `ses_mntd` 返回带类型的结果，包含观测值、零模型均值和标准差、z 分数、秩 p 值、运行次数和零模型。空群落或单物种群落的 MPD/MNTD 观测值为 `NaN`。
- `psv` 返回 `PSVResult`，包含 `psv`、`richness`、`variance` 和 `scale_vcv`；少于两个分类单元的群落的 PSV 和方差为 `NaN`。`psr`、`pse`、`psc` 和 `psd` 返回命名元组，其字段为相应函数文档中的指标名称。
- `expected_pd` 和 `pd_variance` 返回 `ExpectedPDResult`；`richness` 为 `1:tree.ntips`，`expected_pd` 为解析期望值，`variance_pd` 为 `upper_bound` 选择的上界或精确 Picante 方差。
- `pd_calc(tree, tips; variant=...)` 为 tip 子集返回一个标量；矩阵输入返回包含 `pd`、`richness` 和规范化 `variant` 的命名元组。支持 `TBL`、`MST`、`UEH`、`SBL` 和 `TIP`。
- `pd_bootstrap(tree, ntips; ...)` 返回重复计算得到的标量向量；`pd_bootstrap(community, tree; ...)` 返回包含 `observed`、`replicates` 和 `richness` 的 `PDBootstrapResult`。

### 性状结果类型

`pagel_lambda` 返回 `PagelLambdaResult`，包含前述 lambda 拟合值、profile 区间、似然比检验、拟合诊断和模型汇总字段。DataFrame 方法在拟合前按分类单元标签对齐数据。

`phylogenetic_signal_k` 返回 `PhyloSignalResult`，包含 `k`、观测及零模型 PIC 方差、零模型标准差、秩 p 值、z 分数、重排后的性状值、tip 标签、有效重复次数和状态。其零模型打乱 tip 上的性状值；显式传入 `rng` 可复现结果。默认使用 `StableRNG(20260713)`。

`phylogenetic_contrasts` 返回 `PICResult`，包含 `contrasts`、内部节点 ID、预期对比距离、局部祖先值、tip 标签以及 `scaled`/`weighted` 标志。PIC 要求严格二叉树，且枝长有限并严格大于零；多分叉树需显式使用 `multi2di`。PIC 并非通用祖先状态重建 API。

### 边界情况

- `pd` 对空群落给出零 PD；`include_root=true` 时，单物种群落的 PD 为根到 tip 的路径长度。
- 空群落或未定义的群落指标返回 `NaN`，不会静默删除行；UniFrac 和 PhyloSor 使用上述对角线规则。
- 负数、无穷大或 NaN 群落值会触发 `ArgumentError`。
- K/PIC 拒绝缺失和非有限性状值。
- RaoD 拒绝非超度量树。PIC 拒绝多分叉及零或非有限枝长。
- 随机化方法接受外部传入的 `AbstractRNG`；使用种子的 `nri_nti` 在相同种子下保持确定性。

## 性能基准

命令形式：

```julia
tree = read_tree(path)
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

## 功能范围

EvoKits 专注于可移植、以树为中心的工具：

- Newick 读写和 `CompactTree` 转换。
- 树修剪、超度量转换和多分叉解析。
- 树身份标识及 R/ape 映射。
- ED、PD、SES-PD、SES-MPD 和 SES-MNTD 指标。
- Pagel's lambda 和性状对齐。
- 二元性状系统发育信号的 Fritz-Purvis D。
- EDGE2 得分和灭绝概率工具。
- 使用确定性零模型抽样的本地物种池 NRI/NTI。
