# API 约定

[首页](../README_zh.md) · [English](api.md)

组合示例、排序与结果类型约定，以及边界情况。

## 完整示例

以下示例在同一会话中使用树、群落和性状分析功能：

```julia
using EvoKits
using Random

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end

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

### 后验树输入

`read_trees(paths; burnin=0, thin=1)` 返回仅消费一次、逐棵产生 `CompactTree` 的 `TreeStream`；do-block 在提前退出时也关闭文件。`tree_sampler` 返回可迭代的 `TreeSampleResult`，树相互独立，抽样索引、来源路径、文件内记录序号与行号均按抽取顺序排列。筛选分别应用于每个文件，抽样对全部候选记录等概率，默认有放回。格式及报错规则见[树输入说明](trees_zh.md#后验树流式读取与抽样)。

Newick/NEXUS 格式逐文件识别。NEXUS 的树号只统计 `TREE`/`UTREE` 定义，翻译表仅在其所在的 `TREES` 块内生效。

### 共同排序约定

`tree.tip_labels` 是标准 tip 顺序。tip 位置为 `1:tree.ntips` 内的整数；紧凑节点 ID 为 `1:tree.nnodes` 内的整数。接受标签的函数通过 `tree.tipname_to_id` 查找。

群落矩阵以行为群落、列为 tip。所有返回向量或矩阵均保持群落行顺序。`tip_labels` 描述输入列顺序，并将列重排至 `tree.tip_labels` 顺序；它不会重命名或重排群落行。重复、缺失或未知的标签会触发 `ArgumentError`。

群落值必须有限且非负。大于零表示存在；丰度加权方法将正值作为权重。

性状向量默认按树顺序排列。提供 `tip_labels` 时，它描述所提供向量的顺序，程序会将向量重排至 `tree.tip_labels` 顺序。K 和 PIC 拒绝缺失值及非有限值。

### 矩阵语义

`community_longtable(data)` 返回包含 `table`、`sample_labels`、`species_labels` 的命名元组。表的列为 `sample`、`species`、`value`，每个正值一行，先按物种列序、再按群落行序。重建时使用 `CommunityData(table; sample_labels, species_labels)` 并保留两份完整标签列表，以恢复空维度和原顺序。导出保留丰度及数值类型；重建沿用长表规则，包括 Bool 转 Int 计数。参见[往返示例](community_zh.md#导出长表与重建群落)。

`comdist`、`comdistnt`、`unifrac` 和 `phylosor` 返回普通的 `n x n` `Float64` 矩阵，索引顺序与群落行顺序一致。

- `comdist`：两个群落之间的平均成对距离；对角线为包含同种零距离的群落内均值；空群落返回 `NaN`。
- `comdistnt`：默认对两个方向的最近分类单元距离均值等权合并；`abundance_weighted=true` 时按原群落总丰度合并。默认非空对角线为零；`exclude_conspecifics=true` 时为群落内到不同物种的最近距离均值，单物种群落为 `NaN`。排除按物种身份判定，不排除不同物种之间的零距离。
- `unifrac`：非加权枝长距离；对角线为零，包括空群落之间的比较。当非空群落跨越正枝长时，空群落与非空群落之间的距离为一。
- `phylosor`：`2 * shared_PD / (PD1 + PD2)`；非空群落的对角线为一，空群落为零。

`raoD` 返回包含 `Dkk`、`Dkl`、`H`、`total`、`alpha`、`beta` 和 `Fst` 的 `RaoDResult`，使用 `patristic / 2`。树须为超度量树（根到 tip 的绝对容差 1e-8），枝长有限、非负且总和有限。丰度加权默认为 true；关闭时仅将群落内组成改为存在物种等权，合并群落及群落间权重仍按原丰度计算。空群落使 alpha/beta/Fst 为 `NaN`，但合并 total 仍可有效。Dkk、Dkl 和 H 各自拥有独立数组，按输入群落顺序排列，详见 [Rao 分解规则与示例](community_zh.md#rao-多样性分解)。

### 群落结果类型

- `pd` 返回 `PDResult`，包含 `sample_labels`、`pd`、`richness` 和 `include_root`。
- `weighted_endemism` 返回 `WeightedEndemismResult`，包含 `sample_labels`、`we`、`richness` 和按输入物种列顺序排列的 `species_occupancy`。
- `phylo_endemism` 返回 `PhyloEndemismResult`，包含 `sample_labels`、`pe`、`richness`、按树的边顺序排列的 `branch_occupancy` 和 `weighted`。
- `relative_pd` 返回 `RelativePDResult`，包含 `sample_labels`、`rpd`、`pd_normalized`、`pd_equal_normalized`、`richness`、`total_branch_length` 和 `positive_branch_count`。
- `phylobeta` 返回 `PhylobetaResult`，包含 `sample_labels`、`total`、`turnover`、`nestedness`、`index_family` 和 `pairs`。
- `beta_multi` 返回 `BetaMultiResult`，包含标量 `total`/`turnover`/`nestedness`、`index_family`、所选 `sample_labels`、原始行号 `sample_indices`、所选群落 `richness` 和合并物种数 `pooled_richness`，详见[多群落分类 beta](community_zh.md#多群落分类-beta)。
- `phylobeta_multi` 返回 `PhylobetaMultiResult`，包含三个标量成分、`index_family`、所选 `sample_labels`、原始行号 `sample_indices`、`richness`、含根 `pd` 及所选分支并集的 `pooled_pd`；也支持 `CommunityContext`，详见[多群落系统发育 beta](community_zh.md#多群落系统发育-beta)。
- `phylobeta_ses` 返回 `PhylobetaSESResult`，包含群落标签和丰富度、`observed`、`null_mean`、`null_sd`、`z_score`、`p_value`、`valid_replicates`、`pairs`、所选指数族/成分及零模型参数。
- `nmds` 返回 `NMDSResult`，包含标签、`points`、Kruskal `stress`、压缩的输入/坐标/拟合距离、最佳起点历史和终止原因，以及各起点的 stress、迭代数和收敛诊断。
- `align_nmds` 返回 `NMDSAlignmentResult`，包含原群落行序的坐标、共同标签、旋转/尺度/中心、锚点残差及 RMSE。
- `nmds_colors` 返回 `NMDSColorResult`，包含标签、双轴 `points`、对应的 `rgb`/`hex`、可复用的 `NMDSColorScale` 及 `clipped` 标记，详见[坐标对齐与配色](regionalization_zh.md#坐标对齐与配色)。
- `fitgom` 返回 `GoMResult`，包含带标签的 `omega`/`theta`、硬归属、NEF 目标、EM 残差及各起点收敛诊断；`indicators(fit)` 筛选其物种，详见[混合归属](regionalization_zh.md#区域混合归属)。
- `ses_pd`、`ses_mpd` 和 `ses_mntd` 返回带类型的结果，包含观测值、零模型均值和标准差、z 分数、秩 p 值、运行次数和零模型。空群落的 MPD/MNTD 为 `NaN`；单物种群落也为 `NaN`，但正质量的丰度加权 MPD 为零。
- `psv` 返回 `PSVResult`，包含 `psv`、`richness`、`variance` 和 `scale_vcv`。方差对应从全部树 tips 中按该丰富度均匀无放回抽样，完整物种池的方差为零，修正了原公式的叶序依赖。少于两个分类单元的群落的 PSV 和方差为 `NaN`；`compute_var=false` 将方差留为 `NaN`。`psr`、`pse`、`psc` 和 `psd` 返回命名元组，其字段为相应函数文档中的指标名称，详见 [PSV 与 PD 矩](community_zh.md#psvpsr-与解析-pd-矩)。
- `pse` 使用相对丰度，正值等丰度时等于 PSV；`psc` 仅使用存在信息，排除同一物种身份，保留不同物种间的零差异。两者返回相应指标及 `richness`，少于两个物种时为 `NaN`，原始协方差模式的 PSC 可为负。`psd` 按群落顺序返回 `psv`、`var_psv`、`psc`、`psr`、`var_psr`、`pse` 和 `richness`，各字段独立拥有数组，沿用单项指标规则。
- `expected_pd` 和 `pd_variance` 返回 `ExpectedPDResult`；`richness` 为 `1:tree.ntips`。每个 tip 以 `richness/tree.ntips` 的概率独立出现，因此丰富度是期望数量。`expected_pd` 为 PD 均值；`variance_pd` 在 `upper_bound=false` 时为精确 Bernoulli 方差，默认采用 picante 的协方差上界。
- `pd_calc(tree, tips; variant=...)` 为 tip 子集返回一个标量；矩阵输入返回包含 `pd`、`richness` 和规范化 `variant` 的命名元组。支持 `TBL`、`MST`、`UEH`、`SBL` 和 `TIP`。
- `pd_bootstrap(tree, ntips; ...)` 返回重复计算得到的标量向量；`pd_bootstrap(community, tree; ...)` 返回包含 `observed`、`replicates` 和 `richness` 的 `PDBootstrapResult`。

### 分组与保护指标结果

| 函数 | 结果与顺序 | 说明 |
|---|---|---|
| `optimal_functional_groups` | `FunctionalGroupSelectionResult`：递增候选 k、诊断、可为空的建议及可复用分组 | [聚合](aggregation_zh.md) |
| `indicators` | `IndicatorSpeciesResult`：区域 × 排名，原物种索引，候选不足补 missing | [分区](regionalization_zh.md) |
| `phylobuilder` | `PhyloBuilderResult`：独立树；按请求物种顺序的状态与原树锚点 ID | [树操作](trees_zh.md) |
| `edge` | `EDGEResult`：完整树 tip 顺序的 ED/GE/EDGE 及缺失类别状态 | [保护指标](conservation_zh.md) |
| `hotspots`、`coldspots` | `MetricSelectionResult`：原序标签与标记、分位数阈值及有效次数 | [保护指标](conservation_zh.md) |
| `classify_values` | `ValueClassificationResult`：原序类别、上界及可用的 Jenks SSE | [保护指标](conservation_zh.md) |

这些操作均不修改输入。向量指标汇总允许有限负值；`invalid=:skip` 在 missing/NaN 原位置保留 missing 标记或类别。这与群落输入必须非负的约定分别适用。

### 性状结果类型

`pagel_lambda` 返回 `PagelLambdaResult`，包含[系统发育信号](signals_zh.md)中说明的 lambda 拟合值、profile 区间、似然比检验、拟合诊断和模型汇总字段。DataFrame 方法在拟合前按分类单元标签对齐数据。

`phylogenetic_signal_k` 返回 `PhyloSignalResult`，包含 `k`、观测及零模型 PIC 方差、零模型标准差、秩 p 值、z 分数、重排后的性状值、tip 标签、有效重复次数和状态。其零模型打乱 tip 上的性状值；显式传入 `rng` 可复现结果。默认使用 `StableRNG(20260713)`。

`phylogenetic_contrasts` 返回 `PICResult`，包含 `contrasts`、内部节点 ID、预期对比距离、局部祖先值、tip 标签以及 `scaled`/`weighted` 标志。PIC 要求严格二叉树，且枝长有限并严格大于零；多分叉树需显式使用 `multi2di`。PIC 并非通用祖先状态重建 API。

### 边界情况

- `pd` 对空群落给出零 PD；`include_root=true` 时，单物种群落的 PD 为根到 tip 的路径长度。
- 空群落或未定义的群落指标返回 `NaN`，不会静默删除行；UniFrac 和 PhyloSor 使用上述对角线规则。
- 负数、无穷大或 NaN 群落值会触发 `ArgumentError`。
- K/PIC 拒绝缺失和非有限性状值。
- RaoD 拒绝非超度量树。PIC 拒绝多分叉及零或非有限枝长。
- 随机化方法接受外部传入的 `AbstractRNG`；使用种子的 `nri_nti` 在相同种子下保持确定性。
