# 群落多样性与结构

[首页](../README_zh.md) · [English](community.md)

带标签的稠密与稀疏输入、PD、特有性、beta 分解、NRI/NTI、可复用上下文及零模型。

本页示例共用以下小树：

```julia
using EvoKits

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
```

## 群落指标

`CommunityData` 将群落数值、群落标签和物种标签一起保存，接受矩阵、长表 `DataFrame` 或行列索引三元组。矩阵输入会被复制，默认保留数值类型和稠密/CSC 存储；设置 `storage=:sparse` 可转为稀疏存储。长表和三元组直接构建稀疏矩阵，无需先展开成稠密宽表。

```julia
using DataFrames

records = DataFrame(
    sample = ["s1", "s1", "s2", "empty"],
    species = ["A", "C", "A", "A"],
    value = [2, 1, 3, 0],
)
data = CommunityData(records; sample=:sample, species=:species, value=:value)
data.matrix
data.sample_labels   # ["s1", "s2", "empty"]
data.species_labels  # ["A", "C"]

weighted_endemism(data).we
beta_diss(data; pairs=[(1, 2)])
phylo_endemism(data, tree).pe

# 先对齐完整树；需要复用分支集合和占据统计时再准备上下文。
aligned = align_community(data, tree)
ctx = prepare_community(aligned, tree)
pd(ctx).pd
```

长表默认按标签首次出现的顺序排列。可显式提供完整的 `sample_labels`、`species_labels` 控制输出顺序，并保留空群落和未出现物种；列表须包含表中的全部标签。标签接受非空字符串或整数，整数转为字符串，且不能重复。矩阵未提供标签时生成 `sample_1`、`species_1` 等名称，树相关分析需要提供实际物种名。数值必须为有限非负实数，缺失值会报错；设置 `value=nothing` 时按记录次数计数。

重复的群落—物种记录默认求和；设置 `duplicates=:error` 时拒绝重复，包括零值重复记录。求和溢出会报错。输入保留丰度，存在型指标才按正值解释。索引三元组的行列号从一开始，对应显式标签列表：

```julia
data = CommunityData([1, 1, 2], [1, 2, 1], [2, 1, 3];
    sample_labels=["s1", "s2"], species_labels=["A", "C"])
```

树相关调用自动按物种名对齐，为输入中没有的树物种补零列，保留完整树及根。输入中任何树外物种都会报错。`align_community(data, tree; missing_tips=:error)` 还会拒绝缺失的树物种。对齐保留稠密/CSC 存储；先显式对齐一次，后续调用可复用已对齐矩阵，预处理上下文进一步复用分支与占据统计。`CommunityData` 内部数组应视为只读；构造、显式对齐、指标结果和上下文各自拥有独立数组。

现有接受群落矩阵的指标及 `prepare_community` 均支持 `CommunityData`。结果保留原行顺序，已有带标签的结果类型自动使用其群落标签。标签在构造时设置，指标调用中不再另传 `tip_labels` 或 `sample_labels`。SES 的观测与随机化使用稀疏物种集合，也支持丰度加权 MPD/MNTD；其它部分旧指标仍有稠密内部工作区。原有矩阵调用继续可用：

```julia
community = [
    1 0 1 0
    0 1 1 0
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

`ses_pd`、`ses_mpd` 和 `ses_mntd` 只保留一份随机化结果及可复用的物种/树工作区，不再构建稠密“群落 × 物种”随机化缓冲区。支持 `:taxa_labels`、`:richness`、`:frequency`、`:sample_pool`、`:phylogeny_pool`、`:trialswap` 和 `:independentswap`。存储重构保留原抽样规则及基于排列的随机数消耗：richness 保留行和，frequency 保留列和，交换模型同时保留两者。现有 `:sample_pool` 与 `:richness` 一样使用完整树物种池，`:phylogeny_pool` 还会置换物种标签。丰度加权 MPD/MNTD 仍随机化存在状态，并将原丰度固定在树 tip 位置；零权重的随机物种仍作为 MNTD 最近邻候选。

三种 SES 均按有限的随机值计算均值和样本标准差：没有有效随机值时，均值、标准差和 z 均为 `NaN`；只有一个有效随机值时，标准差及 z 为 `NaN`；零模型标准差为零时，z 为 `NaN`。秩 p 值保留原约定 `(1 + 小于观测值的次数 + 0.5 × 相等次数) / (runs + 1)`，非有限随机值不增加秩计数。PD 已修复重复抽样间分支访问标记复用导致的漏计；不包含根时直接累加 MRCA 以下分支，避免用相减去除很长的共有根路径。

SES 虚拟边在深度相减可能丢失精度时，也会回退到原始分支累加，以保留长共享路径下的极短距离。

如需 Caper 兼容的多分叉 ED 处理方式，可使用 `evol_distinct(tree; type = :fair_proportion, polytomy_correction = :isaac)`；其他支持的模式为 `:mooers` 和 `:none`。

## NRI 和 NTI

NRI 和 NTI 描述相对于零模型的群落系统发育结构，分别是 MPD 和 MNTD 标准化效应量的相反数：`NRI = (零模型 MPD 均值 - 观测 MPD) / 零模型 MPD 标准差`，NTI 对应使用 MNTD。正值表示相对于所选零模型的系统发育聚集，负值表示离散。

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

批量分析空间网格时，可传入 `grid_id`，使每个网格的随机种子由全局种子稳定派生，不受处理顺序影响：

```julia
result = nri_nti(tree, pool, observed;
    nreps = 999,
    seed = 20260713,
    grid_id = "1745",
)
```

派生后的实际种子保存在 `result.seed`。不要同时传入 `grid_id` 和显式的 `rng`，因为传入的 RNG 已定义随机数流。

对于大量网格，`nri_nti_table` 接受两张标准分隔长表：一张为本地物种池，一张为观测群落；两表均使用 `grid_id`、`species` 两列。表格可以是 CSV、TSV、TXT 或 CSV.jl 支持的其他分隔格式；需要显式指定制表符时使用 `delim = '\t'`。每个文件中同一网格的记录必须连续。Julia 启动了多个线程时，可并行处理网格。请将下例中的文件路径替换为实际输入。

```julia
batch = nri_nti_table(tree, "local_pool.tsv", "observed_community.tsv";
                        nreps = 999, seed = 20260713, threaded = true,
                        delim = '\t')
batch.grid_ids
batch.results
```

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

## 系统发育 beta SES

`phylobeta_ses` 复用现有 SES 的稀疏随机化和在线统计，对一个系统发育 beta 成分计算标准化效应量。每次选择 `:sorensen` 或 `:jaccard`，以及 `:total`、`:turnover` 或 `:nestedness` 之一。

```julia
using Random

effect = phylobeta_ses(comm, tree; index_family = :sorensen,
    component = :turnover, null_model = :taxa_labels, runs = 999,
    pairs = [(1, 2)], rng = MersenneTwister(42))
effect.observed
effect.null_mean
effect.null_sd
effect.z_score
effect.p_value
effect.valid_replicates
```

支持稠密/CSC 矩阵和 `CommunityData`，所有正丰度视为存在。完整树物种池始终参与抽样，包括未观测物种；`CommunityData` 按物种名对齐并补零，直接矩阵输入须保留每个树 tip 的列。

| 零模型 | 随机化与约束 |
| --- | --- |
| `:taxa_labels`（别名 `:tipshuffle`） | 所有群落共用一次树 tip 标签置换，保留全部分类学群落对重叠数。 |
| `:richness`、`:sample_pool` | 各群落独立从完整树物种池抽取原丰富度的物种。 |
| `:phylogeny_pool` | 丰富度抽样后，再共用一次标签置换。 |
| `:frequency` | 保持每个物种的群落占据数，随机分配其出现位置；群落丰富度可以变化。 |
| `:trialswap`、`:independentswap` | 棋盘格交换同时保持行列和，每次重复从观测群落重新开始。 |

交换次数由 `iterations` 和 `max_attempts_per_swap` 控制，有限交换不保证均匀抽样。对比不同成分时，每次新建相同种子的 RNG 并使用相同参数，即可使用相同随机样本。

`PhylobetaSESResult` 包含上述六组数组、群落标签和丰富度、指定对及零模型参数。默认返回完整对称矩阵；`pairs` 返回按请求顺序排列的向量，保留重复、反向及自配对。每次仍随机化完整输入，仅构建指定端点的分支集合；在线统计不保留各次随机距离矩阵，指定对输出也不分配完整群落间矩阵。

均值和样本 SD 只使用有限随机值，`valid_replicates` 逐对记录次数。没有有效值时均值/SD/z/p 为 `NaN`；只有一个时 SD/z 为 `NaN`；常量零模型保留 SD=0、z=`NaN`。观测 beta 非有限时，即使有有限的零模型统计，z/p 仍为 `NaN`。空群落和零 PD 群落沿用 `phylobeta` 规则，包括未定义的对角线。

下尾中秩分数为 `p = (1 + 小于观测值的次数 + 0.5 × 相等次数) / (有效次数 + 1)`，相等按计算值精确判定。它不是双尾 p 值；分母采用有效次数，旧的群落内 SES 仍保留 `runs + 1`。R 的 `rowwise`、`colwise` 标签重排会被对齐抵消，因此拒绝这两个模型；也不沿用 R 将零 SD 替换为一的做法。SES z 值可能为负或未定义，不宜作为 `phyloregion` 的距离输入。
