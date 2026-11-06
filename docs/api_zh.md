# API 约定

[文档目录](README_zh.md) · [首页](../README_zh.md) · [English](api.md)

本页：[排序](#共同排序约定) · [矩阵](#矩阵语义) · [群落结果](#群落结果类型) · [其他结果](#分组与保护指标结果) · [性状结果](#性状结果类型) · [边界情况](#边界情况)

本页用于查阅排序、结果形状和边界行为。模型与示例见[任务指南](README_zh.md)，算法成本与实测集中在[性能说明](performance_zh.md)。

## 完整示例

[快速开始](../README_zh.md#快速开始)串起带标签群落、PD、PE 和分区。[群落输入](community_input_zh.md)说明构造与对齐，[系统发育信号](signals_zh.md)说明性状分析。

## 公共 API 约定

### 共同排序约定

| 输入或结果 | 顺序 |
| --- | --- |
| 树 | `tree.tip_labels` 为标准 tip 顺序。tip 位置为 `1:tree.ntips`，紧凑节点 ID 为 `1:tree.nnodes`。 |
| 群落矩阵 | 行为群落、列为物种。省略 `tip_labels` 时按树 tip 顺序；提供时描述输入列序，再按树对齐。直接矩阵调用须包含全部树 tips。 |
| `CommunityData` | 按物种名对齐，缺少的树物种补零；未知物种报错，即使对应列全零。群落顺序保持不变。 |
| 性状向量 | 默认按树 tip 顺序；`tip_labels` 可描述另一种输入顺序。 |
| 群落结果 | 向量和矩阵轴保持输入群落顺序；显式选取子集时按请求顺序。 |

群落值须为有限非负实数。正值表示存在；丰度加权方法使用其大小。标签须唯一且非空，整数标识符转为字符串。`tree.tipname_to_id` 用于按标签查找；节点 ID 不等于 tip 位置。

构造和准备函数独立保存数组。上下文字段应视为只读快照，输入改变后需重新准备；返回的指标数组可独立修改。参见[上下文选择](community_input_zh.md#重复分析)和[距离上下文约定](community_distances_zh.md#群落间系统发育距离)。

### 矩阵语义

| 函数 | 默认输出与对角线 | 指定 `pairs` |
| --- | --- | --- |
| `comdist` | `n×n` 跨群落物种对距离均值；对角线包含同种零距离，通常不为零 | 按请求顺序的向量 |
| `comdistnt` | `n×n` 对称最近分类单元均值；默认非空对角线为零，排除同种时为群落内到其他物种的最近距离均值 | 按请求顺序的向量 |
| `unifrac` | `n×n` 分支差异矩阵，对角线为零 | 指定 beta 对可用 `phylobeta` |
| `phylosor` | `n×n` 分支相似度；正 PD 群落对角线为一，零 PD 为零 | 指定 beta 对可用 `phylobeta` |
| `beta_diss`、`phylobeta` | 三份对称 `n×n` 成分数组 | 三份按请求顺序的向量 |
| `phylobeta_ses` | 对称的观测、零模型及 SES 数组 | 按请求顺序的向量 |
| `mean_dist` | 对称的距离及有限观测次数 | 按请求顺序的向量 |

`pairs` 使用完整输入群落行号，保留重复、反向和自身对。空请求返回空向量，仍校验输入。矩阵型数值输出为普通 `Float64` 数组，计数数组为整数。定义见[群落距离](community_distances_zh.md)、[beta 分解](community_zh.md#特有性相对-pd-与-beta-分解)及[零模型](null_models_zh.md)。

`community_longtable(data)` 返回 `(table, sample_labels, species_labels)`；表含 `sample`、`species`、`value` 三列，每个正单元格一行，先按物种列序、再按群落行序。保留两份完整标签列表，才能在[重建](community_input_zh.md#导出长表与重建群落)时恢复空群落和全零物种。

### 群落结果类型

| 函数 | 结果 | 主要字段 |
| --- | --- | --- |
| `pd` | `PDResult` | `sample_labels`、`pd`、`richness`、`include_root` |
| `weighted_endemism` | `WeightedEndemismResult` | `sample_labels`、`we`、`richness`、`species_occupancy`（输入物种序） |
| `phylo_endemism` | `PhyloEndemismResult` | `sample_labels`、`pe`、`richness`、`branch_occupancy`（树边序）、`weighted` |
| `relative_pd` | `RelativePDResult` | `sample_labels`、`rpd`、`pd_normalized`、`pd_equal_normalized`、`richness`、`total_branch_length`、`positive_branch_count` |
| `psv` | `PSVResult` | `psv`、`richness`、`variance`、`scale_vcv` |
| `psr`、`pse`、`psc` | 命名元组 | 相应指标数组及 `richness`；PSR 还包含 `variance` |
| `psd` | 命名元组 | `psv`、`var_psv`、`psc`、`psr`、`var_psr`、`pse`、`richness` |
| `expected_pd`、`pd_variance` | `ExpectedPDResult` | `richness`、`expected_pd`、`variance_pd`、`upper_bound` |
| `raoD` | `RaoDResult` | 群落内 `Dkk`、群落对 `Dkl`/`H`、合并 `total`、`alpha`、`beta`、`Fst` |
| `pd_calc` | 树/tip 子集返回标量；群落返回命名元组 | `pd`、`richness`、规范化 `variant` |
| `pd_bootstrap` | 树/丰富度返回向量；群落返回 `PDBootstrapResult` | `observed`、群落 × 重复次数的 `replicates`、`richness` |
| `beta_diss`、`phylobeta` | `BetaDissResult`、`PhylobetaResult` | `sample_labels`、`total`、`turnover`、`nestedness`、`index_family`、`pairs` |
| `beta_multi` | `BetaMultiResult` | 标量成分、`index_family`、所选 `sample_labels`、原行号 `sample_indices`、`richness`、`pooled_richness` |
| `phylobeta_multi` | `PhylobetaMultiResult` | 标量成分、`index_family`、所选 `sample_labels`、原行号 `sample_indices`、`richness`、含根 `pd`、所选分支并集的 `pooled_pd` |
| `ses_pd`、`ses_mpd`、`ses_mntd` | 带类型的群落 SES 结果 | `sample_labels`、`ntaxa`、带前缀的观测/零模型/秩/z/p 数组及运行设置；MPD/MNTD 另含 `mpd_nri`/`mntd_nti` |
| `phylobeta_ses` | `PhylobetaSESResult` | `sample_labels`、`richness`、`observed`、`null_mean`、`null_sd`、`z_score`、`p_value`、`valid_replicates`、`pairs` 及运行设置 |

`expected_pd` 默认返回 `1:tree.ntips`；`richness` 接受该范围内的整数、整数向量或区间。三个结果数组均保留请求顺序和重复值，包括单元素和空请求。其独立出现模型与 PSV 方差的固定样本量模型不同，详见[多样性矩](community_zh.md#psvpsr-与解析-pd-矩)。

`raoD(...; pairwise=false)` 返回独立的 `0×0` `Dkl`/`H`；`Dkk` 和四个标量定义不变，详见 [Rao 分解规则](community_zh.md#rao-多样性分解)。多群落 beta 是整体标量分解，不是两两 beta 的均值。

### 后验树输入

`read_trees(paths; burnin=0, thin=1)` 返回逐棵产生 `CompactTree`、仅消费一次的 `TreeStream`。使用 do-block 可在提前退出时关闭文件。筛选及 Newick/NEXUS 识别均逐文件进行；NEXUS 树号只统计 `TREE`/`UTREE` 定义，翻译表仅在其 `TREES` 块内生效。

`tree_sampler` 返回可迭代的 `TreeSampleResult`，包含独立的 `trees` 及按抽取顺序排列的 `sample_indices`、`source_paths`、`source_indices`、`source_lines`，另含 `population_size`、`replace`、`burnin`、`thin`。对全部候选记录等概率抽样，默认有放回。见[后验树输入](trees_zh.md#后验树流式读取与抽样)。

### 后验距离均值

`mean_dist(data, trees; index_family=:sorensen, component=:total, pairs=nothing, nan_policy=:skip, block_size=256)` 直接消费树迭代器。矩阵未传 `tip_labels` 时，由第一棵树确定列含义；`CommunityData` 自带物种标签。单参数形式汇总已有距离或 beta 结果，并按群落标签对齐。

`MeanDistResult` 包含 `sample_labels`、`distances`、`counts`、`pairs`、`replicates`、`component`、`index_family`。`counts` 按群落对统计有限值，`replicates` 统计全部输入；无有效值返回 NaN。`nan_policy=:propagate` 在任一重复未定义时保留 NaN。见[后验均值](regionalization_zh.md#多棵后验树的平均距离)。

### 分组与保护指标结果

| 函数 | 结果与顺序 | 指南 |
| --- | --- | --- |
| `optimal_functional_groups` | `FunctionalGroupSelectionResult`：递增候选 k、诊断、可为空的建议及可复用分组 | [聚合](aggregation_zh.md) |
| `nmds` | `NMDSResult`：原群落序 `points`、Kruskal `stress`、压缩输入/坐标/拟合距离、最佳起点和各起点诊断 | [NMDS](regionalization_zh.md#非度量多维尺度分析nmds) |
| `align_nmds` | `NMDSAlignmentResult`：原序坐标、共同标签、旋转/尺度/中心、锚点残差及 RMSE | [对齐](regionalization_zh.md#坐标对齐与配色) |
| `nmds_colors` | `NMDSColorResult`：标签、双轴 `points`、`rgb`/`hex`、可复用 `NMDSColorScale`、`clipped` | [配色](regionalization_zh.md#坐标对齐与配色) |
| `fitgom` | `GoMResult`：带标签的群落 × k `omega`、物种 × k `theta`、硬归属、目标值及收敛诊断 | [混合归属](regionalization_zh.md#区域混合归属) |
| `indicators` | `IndicatorSpeciesResult`：区域 × 排名，原物种索引，候选不足补 missing | [指示物种](regionalization_zh.md#区域指示物种) |
| `phylobuilder` | `PhyloBuilderResult`：独立树，按请求物种顺序的状态和原树锚点 ID | [树操作](trees_zh.md) |
| `edge` | `EDGEResult`：完整树 tip 序、ED/GE/EDGE 及缺失类别状态 | [保护指标](conservation_zh.md) |
| `hotspots`、`coldspots` | `MetricSelectionResult`：原序标签/标记、分位数阈值及有效次数 | [保护指标](conservation_zh.md) |
| `classify_values` | `ValueClassificationResult`：原序类别、上界及可用的 Jenks SSE | [保护指标](conservation_zh.md) |

GoM 权重描述组成，不代表不确定性。GoM 或 NMDS 收敛不等于达到全局最优。向量指标汇总允许有限负值，`invalid=:skip` 在 missing/NaN 位置保留 missing 标记/类别；群落丰度仍须非负。

### 性状结果类型

| 函数 | 结果及主要字段 |
| --- | --- |
| `pagel_lambda` | `PagelLambdaResult`：lambda 拟合值、profile 区间、似然比检验、拟合诊断和模型汇总；DataFrame 输入按分类单元标签对齐 |
| `phylogenetic_signal_k` | `PhyloSignalResult`：`k`、观测/零模型 PIC 方差、零模型 SD、秩 p 值、z 分数、对齐性状/tip 标签、有效次数和状态 |
| `phylogenetic_contrasts` | `PICResult`：对比值、内部节点 ID、预期距离、局部祖先值、tip 标签、`scaled`/`weighted` |

K 打乱 tip 性状，默认 RNG 为 `StableRNG(20260713)`。PIC 要求严格二叉树及有限正枝长，多分叉需显式使用 `multi2di`；PIC 并非通用祖先状态重建 API。见[性状模型及诊断](signals_zh.md)。

### 边界情况

| 情形 | 行为 |
| --- | --- |
| 空群落 / 单物种 PD | `pd`：空群落为零；单物种在 `include_root=true` 时为根到 tip 长度，false 时为 NaN。`pd_calc` 各变体按分支规则计算，单物种 MST 为零。 |
| 空群落 MPD/MNTD | NaN；单物种也为 NaN，但有正质量的丰度加权 MPD 为零 |
| 少于两个物种 | PSV/PSR/PSE/PSC 及 PSV/PSR 方差为 NaN |
| 空群落/零 PD 的 beta | `beta_diss` 按丰富度、`phylobeta` 按 PD 判断。两者均为零时，全部成分为 NaN；仅一者为零时，总差异为一，周转/嵌套为 NaN。 |
| UniFrac/PhyloSor 分母为零 | 均返回零，包括非空但 PD 为零的群落 |
| 非法群落值 | 负数、无穷大或 NaN 触发 `ArgumentError`，不静默删除行 |
| 树限制 | RaoD 要求超度量树，根到 tip 绝对容差为 `1e-8`；相关矩阵模式的 PSV 系列要求根到 tip 长度为正 |
| 缺失/非有限性状 | K/PIC 拒绝输入 |
| 随机化方法 | 显式传入 `AbstractRNG` 以复现；`nri_nti` 也接受稳定种子。不同函数的默认 RNG 不完全相同。 |

SES 零 SD 和下尾 p 值分母见[零模型](null_models_zh.md)。聚类/NMDS 要求完整、有限、非负、对称且对角线为零的距离，见[分区](regionalization_zh.md)。PhyloSor 为相似度，comdist 的对角线通常非零，SES z 分数也不能直接作为差异输入。
