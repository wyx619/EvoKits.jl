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
- `weighted_endemism` 返回 `WeightedEndemismResult`，包含 `sample_labels`、`we`、`richness` 和按输入物种列顺序排列的 `species_occupancy`。
- `phylo_endemism` 返回 `PhyloEndemismResult`，包含 `sample_labels`、`pe`、`richness`、按树的边顺序排列的 `branch_occupancy` 和 `weighted`。
- `relative_pd` 返回 `RelativePDResult`，包含 `sample_labels`、`rpd`、`pd_normalized`、`pd_equal_normalized`、`richness`、`total_branch_length` 和 `positive_branch_count`。
- `phylobeta` 返回 `PhylobetaResult`，包含 `sample_labels`、`total`、`turnover`、`nestedness`、`index_family` 和 `pairs`。
- `phylobeta_ses` 返回 `PhylobetaSESResult`，包含群落标签和丰富度、`observed`、`null_mean`、`null_sd`、`z_score`、`p_value`、`valid_replicates`、`pairs`、所选指数族/成分及零模型参数。
- `nmds` 返回 `NMDSResult`，包含标签、`points`、Kruskal `stress`、压缩的输入/坐标/拟合距离、最佳起点历史和终止原因，以及各起点的 stress、迭代数和收敛诊断。
- `ses_pd`、`ses_mpd` 和 `ses_mntd` 返回带类型的结果，包含观测值、零模型均值和标准差、z 分数、秩 p 值、运行次数和零模型。空群落的 MPD/MNTD 为 `NaN`；单物种群落也为 `NaN`，但正质量的丰度加权 MPD 为零。
- `psv` 返回 `PSVResult`，包含 `psv`、`richness`、`variance` 和 `scale_vcv`；少于两个分类单元的群落的 PSV 和方差为 `NaN`。`psr`、`pse`、`psc` 和 `psd` 返回命名元组，其字段为相应函数文档中的指标名称。
- `expected_pd` 和 `pd_variance` 返回 `ExpectedPDResult`；`richness` 为 `1:tree.ntips`，`expected_pd` 为解析期望值，`variance_pd` 为 `upper_bound` 选择的上界或精确 Picante 方差。
- `pd_calc(tree, tips; variant=...)` 为 tip 子集返回一个标量；矩阵输入返回包含 `pd`、`richness` 和规范化 `variant` 的命名元组。支持 `TBL`、`MST`、`UEH`、`SBL` 和 `TIP`。
- `pd_bootstrap(tree, ntips; ...)` 返回重复计算得到的标量向量；`pd_bootstrap(community, tree; ...)` 返回包含 `observed`、`replicates` 和 `richness` 的 `PDBootstrapResult`。

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
