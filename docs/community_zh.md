# 群落多样性与结构

[文档目录](README_zh.md) · [首页](../README_zh.md) · [English](community.md)

本页：[PD 与 bootstrap](#pd-变体与-bootstrap) · [多样性矩](#psvpsr-与解析-pd-矩) · [Rao](#rao-多样性分解) · [特有性与两两 beta](#特有性相对-pd-与-beta-分解) · [多群落分类 beta](#多群落分类-beta) · [多群落系统发育 beta](#多群落系统发育-beta)

按分析问题选择多样性或 beta 指标，计算并解释结果。本页示例共用以下树和带标签群落：

```julia
using EvoKits, SparseArrays

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
community = [1 0 1 0; 0 1 1 0]
data = CommunityData(sparse(community);
    sample_labels = ["west", "east"], species_labels = tree.tip_labels)
```

## 群落指标

| 分析问题 | 函数 | 说明 |
| --- | --- | --- |
| 群落包含多少分支历史？ | `pd`、`pd_calc`、`pd_bootstrap` | [PD 变体](#pd-变体与-bootstrap) |
| 物种有多大变异、均匀性或聚集程度？ | `psv`、`psr`、`pse`、`psc`、`psd` | [多样性矩](#psvpsr-与解析-pd-矩) |
| 独立出现模型下，PD 的期望值和方差是多少？ | `expected_pd`、`pd_variance` | [解析 PD 矩](#psvpsr-与解析-pd-矩) |
| 如何分解 Rao 多样性？ | `raoD` | [Rao 分解](#rao-多样性分解) |
| 哪些群落包含分布局限的物种或分支？ | `weighted_endemism`、`phylo_endemism`、`relative_pd` | [特有性与两两 beta](#特有性相对-pd-与-beta-分解) |
| 两个群落有何差异？ | `beta_diss`、`phylobeta` | [两两 beta](#特有性相对-pd-与-beta-分解) |
| 整个群落集合的多样性有多大？ | `beta_multi`、`phylobeta_multi` | [多群落 beta](#多群落分类-beta) |

### 导出长表与重建群落

构造方法、稀疏存储和完整长表往返示例见[群落输入](community_input_zh.md#导出长表与重建群落)。

### 树对齐与指标调用

缺失物种、标签顺序及准备函数选择见[对齐与上下文](community_input_zh.md#树对齐与指标调用)。

### PD 变体与 bootstrap

`pd_calc(community, tree; variant=:tbl)` 按群落顺序返回 `pd`、`richness` 和规范化的变体名；`pd_calc(tree, tips; variant=...)` 为不重复的 tip 名称或位置返回标量。稠密/CSC 矩阵及 `CommunityData` 均按正丰度解释存在。

| 变体 | 纳入的原始分支 |
| --- | --- |
| `:tbl` | 观测 tips 到原树根的全部路径。 |
| `:mst` | 观测 MRCA 下方的连接分支，排除其上方共有根枝。 |
| `:ueh` | 恰有一个观测后代的分支。 |
| `:sbl` | 至少有两个观测后代的分支。 |
| `:tip` | 观测物种在原树上的末端分支。 |

空群落五项均为零；单物种 MST 为零，UEH 包含整条到根路径，TIP 只包含原始末端边。

```julia
using Random

variants = pd_calc(data, tree; variant = :ueh)
draws = pd_bootstrap(data, tree; reps = 99, rng = MersenneTwister(42))
draws.observed
draws.replicates              # 群落 × 重复次数
fixed_richness = pd_bootstrap(tree, 2; reps = 99, method = :MST,
    tip_weights = [1.0, 2.0, 3.0, 4.0], rng = MersenneTwister(42))
```

bootstrap 从完整树 tip 池加权无放回抽样，保留各群落丰富度。权重按树 tip 顺序，转换后须为有限非负 Float64 且至少一项为正，省略时均匀抽样。树输入返回向量；群落输入返回 `PDBootstrapResult`，空群落为零，`include_root=false` 时单物种群落为 NaN，与 `pd` 一致。

若请求丰富度超过正权重 tip 数，剩余位置按原 tip 顺序从零权重物种中填入，与 R 拒绝支持数不足的行为不同。各权重有限但求和溢出时，先缩放再归一化。正权重支持数充足时，R 与 Julia 对齐加权抽样分布，不要求跨语言随机流相同。耗时与存储见[性能说明](performance_zh.md#pd-变体与-bootstrap)。

### PSV、PSR 与解析 PD 矩

`psv` 是物种对差异 `(Cii+Cjj)/2-Cij` 的均值。C 默认为 Brownian 相关矩阵；`scale_vcv=false` 使用协方差，此时该差异等于 `patristic / 2`。`psr` 为 PSV 乘群落丰富度。两者直接使用树上统计量，方差也无需 VCV，接受稠密/CSC 输入和 `CommunityData`。

PSV 方差按完整树物种池在指定丰富度下**均匀无放回**抽样定义，物种池包含未观测 tip；抽取完整物种池时方差为零。与参考 picante 实现中的叶序依赖公式不同，这里的方差不随 tip 排序改变，原始协方差模式也计入不同 tip 方差的贡献。空群落和单物种群落返回 `NaN`，`compute_var=false` 时方差为 `NaN`；相关矩阵模式要求根到每个 tip 的长度为正。

```julia
moment_data = CommunityData(sparse([1 1 0 0; 1 1 1 1]);
    species_labels=["A", "B", "C", "D"])
variability = psv(moment_data, tree)
variability.psv       # [1/2, 5/6]
variability.variance  # [1/18, 0]
richness_metric = psr(moment_data, tree)
moments = expected_pd(tree; upper_bound=false)
moments.expected_pd  # [1.875, 3.5, 4.875, 6.0]
moments.variance_pd  # [2.3671875, 2.375, 1.2421875, 0.0]
selected = expected_pd(tree; richness=[3, 1, 4], upper_bound=false)
selected.richness     # [3, 1, 4]
selected.expected_pd  # [4.875, 1.875, 6.0]
diversity = psd(moment_data, tree)
diversity.pse         # [1/2, 5/6]：等丰度时 PSE = PSV
diversity.psc         # [1/2, 1/2]
pse(sparse([2.0 1 0 0]), tree).pse  # [4/9]
```

`pse` 按相对丰度加权物种对差异；`psc` 为各存在物种到其它物种的最小 `1-Cij` 的均值。原始协方差模式下 PSC 是一减最大共有枝长的均值，可为负，且依赖枝长单位；相关矩阵模式在非超度量树上选中的邻居可能与 patristic 最近物种不同。不同物种间的零长度路径保留。`psd` 一次返回四项指标、PSV/PSR 方差和丰富度，沿用单项调用的定义及 NaN 规则，接受并保留稠密/CSC 与 `CommunityData` 输入；算法及分配开销见[独立性能说明](performance_zh.md#解析多样性矩)。

PSV/PSR（含 PSD 对应输出）的方差始终基于完整原树物种池，包括未观测物种；算法与存储见[性能说明](performance_zh.md#解析多样性矩)。

`expected_pd` 及别名 `pd_variance` 保留另一种模型：每个 tip 以 `richness/tree.ntips` 的概率**独立出现**，richness 是期望物种数。`upper_bound=false` 给出精确 Bernoulli 方差，默认 true 给出 picante 的协方差上界；并非固定样本量稀释曲线。

省略 `richness` 返回完整 `1:tree.ntips` 曲线；指定该区间内的整数、整数向量或区间时，仅计算请求点，保留顺序与重复值。标量仍返回单元素数组，空向量返回空数组；Bool、非整数类型及越界值报错。结果独立拥有数组。即使请求为空，也要求枝长有限、非负且总和有限。[算法与扩展性](performance_zh.md#解析多样性矩)单独说明。

### Rao 多样性分解

`raoD` 要求超度量树，使用 `patristic / 2`，遵循 Picante 的分解约定。返回的 `RaoDResult` 包含群落内 `Dkk`、群落间 `Dkl`、`H[i,j] = Dkl[i,j] - (Dkk[i]+Dkk[j])/2`、合并群落的 `total`、按群落权重加权的 `alpha`、`beta = total-alpha` 和 `Fst = beta/total`。

```julia
rao_data = CommunityData(sparse([2.0 1 0; 0 1 3]);
    sample_labels=["west", "east"], species_labels=["A", "B", "C"])
rao = raoD(rao_data, tree)
rao.Dkk           # [4/9, 3/4]
rao.total         # 8/7
rao.Fst           # 11/24
raoD(rao_data, tree; abundance_weighted=false)
rao_summary = raoD(rao_data, tree; pairwise=false)
rao_summary.Fst    # 11/24
size(rao_summary.Dkl), size(rao_summary.H) # ((0, 0), (0, 0))
```

默认 `pairwise=true` 返回两份群落对矩阵。只需要群落内多样性和合并分解时，设置 `pairwise=false` 跳过群落对计算，仍返回 `RaoDResult`，其中 `Dkl`、`H` 为各自独立的 `0×0` 矩阵，`Dkk`、`total`、`alpha`、`beta`、`Fst` 沿用原定义。

默认按丰度加权。`abundance_weighted=false` 只将各群落内部的存在物种改为等权；合并群落的组成和群落间权重仍使用原丰度，因此 beta/Fst 可能为负。这是现有 Julia 选项，Picante 没有对应参数。

空群落的 Dkk 及 Dkl/H 整行整列为 `NaN`；即使群落权重为零，也使 alpha/beta/Fst 为 `NaN`，但 total 仍可有效。全空输入的 total 为 `NaN`；零输入行返回空数组和 alpha=0。总多样性为零时 Fst 未定义，beta 也为零时返回 `NaN`。稠密/CSC 矩阵和 `CommunityData` 均遵循这些规则，输入不被修改；[算法与存储开销](performance_zh.md#rao-多样性分解)单独说明。

## 特有性、相对 PD 与 beta 分解

这些函数接受普通稠密矩阵和 `SparseMatrixCSC` 稀疏矩阵。行表示群落，列表示物种，数值必须有限且非负，所有正值均按存在处理。直接传矩阵给使用树的函数时，要求每个 tip 对应一列；未出现在任何群落中的树物种也保留为零列。输入列顺序与 `tree.tip_labels` 不同时，通过 `tip_labels` 指明列标签。使用 `CommunityData` 时自动完成物种对齐及补零。

```julia
using SparseArrays

comm = sparse(community)
we = weighted_endemism(comm)
pe = phylo_endemism(comm, tree; tip_labels = tree.tip_labels)
strict_pe = phylo_endemism(comm, tree; weighted = false)
rpd = relative_pd(comm, tree)

we.we
pe.pe
strict_pe.pe
rpd.rpd
rpd.pd_normalized
rpd.pd_equal_normalized

# 分类 beta 使用物种计数，不需要树。
tax_beta = beta_diss(comm; index_family = :sorensen)
tax_beta.total
tax_beta.turnover
tax_beta.nestedness
tax_selected = beta_diss(comm; index_family = :jaccard, pairs = [(1, 2)])

beta = phylobeta(comm, tree; index_family = :sorensen)
beta.total
beta.turnover
beta.nestedness

# 只计算指定群落对，避免分配三个完整的群落间矩阵。
selected = phylobeta(comm, tree; index_family = :jaccard, pairs = [(1, 2)])
selected.total
```

加权特有性 WE 对群落内物种的占据群落数取倒数后求和。加权系统发育特有性 PE 对群落包含的各条分支，按“枝长 / 包含该分支的群落数”求和；设置 `weighted=false` 时计算严格 PE，只累加仅在一个群落出现的分支长度。占据范围基于完整输入，分别对几个批次独立计算 WE 或 PE 会改变统计量。空群落返回零。

相对系统发育多样性 RPD 使用拓扑相同、所有正枝长相等且保留零枝长的参照树。原树和参照树分别归一化至总枝长为一，再计算 PD 比值。本节使用树的指标都包含至输入树根的路径，计算不会修改原树。RPD 的分母为零时返回 `NaN`。

`beta_diss` 使用共享和独有物种数计算分类 beta 分解，返回 `BetaDissResult`，不需要树，也不按丰度大小加权；`phylobeta` 则使用共享和独有枝长。两者均支持 Sørensen 和 Jaccard，返回总差异 `total`、周转 `turnover` 和嵌套 `nestedness`。默认字段为对称矩阵；提供 `pairs` 后返回与指定群落对顺序一致的向量，保留重复、反向和自配对。定义有效时，总差异等于周转与嵌套之和。

`beta_diss` 比较空群落与非空群落时，总差异为一，周转和嵌套为 `NaN`；两个空群落的三项均为 `NaN`，包括对角线，非空群落的对角线为零。`phylobeta` 比较零 PD 群落与正 PD 群落时，总差异为一，周转和嵌套未定义，返回 `NaN`；两个零 PD 群落的三项均为 `NaN`。正 PD 群落的对角线为零。

`weighted_endemism`、`phylo_endemism`、`relative_pd`、`beta_diss` 和 `phylobeta` 都接受 `sample_labels` 并保留群落行顺序，结果在适用时提供丰富度或占据数。`pd`、`unifrac` 和 `phylosor` 也支持稀疏输入处理，不会强制构造稠密群落矩阵。WE、PE 和 PD 无需保存所有 tip 到根的完整路径。树相关指标要求枝长有限且非负，总枝长也必须有限。

对同一份树和群落数据重复计算时，可先准备上下文：

```julia
ctx = prepare_community(comm, tree; tip_labels = tree.tip_labels)
pd(ctx).pd
pd(ctx; include_root = false).pd
weighted_endemism(ctx).we
phylo_endemism(ctx).pe
phylo_endemism(ctx; weighted = false).pe
relative_pd(ctx).rpd
unifrac(ctx)
phylosor(ctx)

# 每次请求都复用标签匹配和分支集合。
first_batch = phylobeta(ctx; pairs = [(1, 2)])
second_batch = phylobeta(ctx; index_family = :jaccard, pairs = [(2, 1)])
```

`prepare_community` 返回 `CommunityContext`，独立保存标签匹配后的分支集合、完整输入中的物种和分支占据统计，以及各群落的汇总值。后续群落对索引始终对应准备时的输入行；WE/PE 的占据范围始终是准备时的完整集合，物种占据数按原输入列顺序返回。群落标签通过准备时的 `sample_labels` 设置，也支持 `pd_dataframe(ctx)`。

原矩阵、树或标签后续发生变化时，上下文仍使用准备时的快照；分析新数据需要重新准备。可以独立修改返回结果；上下文内部数组应视为只读。准备过程会保留分支集合及汇总值，因此有一次性的时间和内存成本，适合同一数据的重复分析。上下文不保留原始树和矩阵，也不缓存完整群落间矩阵；完整矩阵调用仍执行两两计算并分配输出，指定群落对调用仅分配请求的结果。

## 多群落分类 beta

`beta_multi` 按多群落 Sørensen 或 Jaccard 定义，返回所选群落整体的标量 `total`、`turnover` 和 `nestedness`。

```julia
collection = CommunityData(sparse([1 1 1 0; 1 0 0 1; 1 0 0 0]);
    sample_labels = ["s1", "s2", "s3"], species_labels = ["A", "B", "C", "D"])
overall = beta_multi(collection)
overall.total             # 3/5
overall.turnover          # 1/3
overall.nestedness        # 4/15
overall.richness          # [3, 2, 1]
overall.pooled_richness   # 4

subset = beta_multi(collection; samples = ["s3", "s1"], index_family = :jaccard)
subset.sample_labels      # ["s3", "s1"]
subset.sample_indices     # [3, 1]，对应原始输入行
subset.pooled_richness    # 3，基于所选群落重新计数
```

支持稠密/CSC 矩阵和 `CommunityData`，所有正值均按存在处理，丰度大小不参与加权。矩阵标签通过 `sample_labels` 设置，须唯一且非空，默认生成 `sample_1` 等标签。`samples` 接受由不同整数行号或字符串标签组成的向量、元组，并保留请求顺序；整数标签会转成字符串，按标签选择时应传字符串。不接受 Bool 掩码、重复群落或行号与标签混合的选择。输入数值整体校验，包括未选行。结果独立保存标签、原始行号及所选群落的丰富度。

令 A 为所选群落丰富度之和减去合并物种数，M/X 为全部无序群落对独有物种数的最小值/最大值之和。Sørensen 总差异为 `(M+X)/(2A+M+X)`，替代为 `M/(A+M)`；Jaccard 总差异为 `(M+X)/(A+M+X)`，替代为 `2M/(A+2M)`。嵌套采用等价的乘积公式保留很小的数值。这些是整体指标，对 `beta_diss` 的元素求均值得到的是另一个统计量。

至少选择两个群落。空群落保留，全部为空或零物种列时三项均为 `NaN`。只有一个非空群落时，总差异为一，替代与嵌套为 `NaN`；无共享物种且至少两个非空群落时，总差异和替代为一、嵌套为零。其余情况按各成分的分母判断。输入不修改，也不构造完整两两距离矩阵；算法说明见[性能专题](performance_zh.md#多群落分类-beta)。

## 多群落系统发育 beta

`phylobeta_multi` 将相同的多群落分解扩展到分支长度，支持稠密/CSC 矩阵和树、`CommunityData` 和树，或已准备的 `CommunityContext`。沿用上面的 `collection`：

```julia
phylogeny = mktempdir() do dir
    path = joinpath(dir, "phylogeny.nwk")
    write_tree(path, "((A:1,B:2):3,(C:4,D:5):6);")
    read_tree(path)
end
phylo_overall = phylobeta_multi(collection, phylogeny)
phylo_overall.total          # 17/31
phylo_overall.turnover       # 5/19
phylo_overall.nestedness     # 168/589
phylo_overall.pd             # [16, 15, 4]
phylo_overall.pooled_pd      # 21

phylo_context = prepare_community(collection, phylogeny)
phylo_subset = phylobeta_multi(phylo_context;
    samples = ["s3", "s1"], index_family = :jaccard)
phylo_subset.sample_indices  # [3, 1]
phylo_subset.pd              # [4, 16]
phylo_subset.pooled_pd       # 16，仅包含所选群落占据的分支
```

`samples` 沿用 `beta_multi` 的选择和校验规则。`PhylobetaMultiResult` 独立保存所选标签、原始行号、丰富度和含根 PD，三个 beta 成分为标量。正丰度表示存在；矩阵列按树 tip 顺序或 `tip_labels` 匹配，`CommunityData` 按物种标签对齐并为缺失 tip 补零。

上述公式中的 A 改为 PD 之和减去合并 PD，M/X 使用独有分支长度。合并 PD 是所选群落实际占据的分支并集，保留到原树根的路径，未观测物种不增加分支；这与 betapart 的直接树计算一致。本地 phyloregion 核心将整棵树长作为合并 PD，在存在未观测分支时可能得到负的多群落成分，此处不沿用该行为。

上下文调用复用分支集合，并在所选范围内重新统计占据数，不改变 WE/PE 使用的完整群落占据统计。至少选择两个群落；空群落和零 PD 群落均保留。所有 PD 为零时三项均为 `NaN`；只有一个正 PD 群落时，总差异为一，周转和嵌套为 `NaN`；其余情况按各成分分母判断。选择两个群落时退化为 `phylobeta`。输入不修改，也不建立完整两两矩阵；算法说明见[性能专题](performance_zh.md#多群落系统发育-beta)。


### 群落间系统发育距离

`comdist`、`comdistnt`、UniFrac/PhyloSor、指定群落对和重复距离计算见[群落距离](community_distances_zh.md)。

## NRI 和 NTI

本地物种池、正负号约定和批量 NRI/NTI 见[零模型与群落结构](null_models_zh.md#nri-和-nti)。

## 系统发育 beta SES

随机化模型、有效重复数及 z/p 值解释见 [beta SES](null_models_zh.md#系统发育-beta-ses)。
