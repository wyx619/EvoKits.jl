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

### 导出长表与重建群落

`community_longtable(data)` 返回命名元组，包含记录表 `table`（`DataFrame`）和完整的 `sample_labels`、`species_labels` 向量。每个正值群落—物种组合对应一行，列名固定为 `sample`、`species`、`value`；先按物种列序，再按群落行序排列。丰度与小数权重保留原值和数值类型，整数丰度保存在数值列中，零值不导出。

```julia
using SparseArrays

export_source = CommunityData(sparse([0.0 2.5 0; 0 0 0; 3 0.25 0]);
    sample_labels=["s2", "empty", "s1"], species_labels=["C", "A", "absent"])
records = community_longtable(export_source)
records.table
restored = CommunityData(records.table;
    sample_labels=records.sample_labels, species_labels=records.species_labels)
restored.matrix == export_source.matrix # true
restored.sample_labels                 # ["s2", "empty", "s1"]
restored.species_labels                # ["C", "A", "absent"]
```

保存记录表时须同时保存两份标签列表：只有非零记录无法保留空群落、全零物种列，也可能改变群落顺序。全零矩阵、没有行或列的矩阵同样适用。输出数组与输入独立；重建采用 CSC 存储，沿用长表的数值规则，包括将 Bool 转为整数计数。CSC 导出直接遍历存储项，[存储开销](performance_zh.md#群落长表导出)单独说明。

### 树对齐与指标调用

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

矩阵的行表示群落，列表示树的 tip。如果输入列的顺序与树 tip 顺序不同，请提供 `tip_labels`。当距离或多样性未定义时，空群落返回 `NaN`；PD、UniFrac、PhyloSor 和 RaoD 接受单物种群落。

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

bootstrap 使用指数随机键从完整树 tip 池加权无放回抽样，保留各群落丰富度。权重按树 tip 顺序，转换后须为有限非负 Float64 且至少一项为正，省略时均匀抽样。树输入返回向量；群落输入返回 `PDBootstrapResult`，空群落为零，`include_root=false` 时单物种群落为 NaN，与 `pd` 一致。

权重总和有限时，选中子集及 RNG 消耗保持原完整排序实现的结果；并列按原 tip 位置，零权重随机键为 `-Inf`。若丰富度超过正权重支持数，保留旧实现按索引填入剩余 `-Inf` 项的规则，与 R 拒绝支持数不足的行为不同。各权重有限但求和溢出时，先缩放再归一化。R 与 Julia 对齐抽样分布，不要求跨语言随机流相同。稀疏工作区和有限候选堆替代稠密随机群落及完整排序，耗时和剩余成本见[性能说明](performance_zh.md#pd-变体与-bootstrap)。

### PSV、PSR 与解析 PD 矩

`psv` 是物种对差异 `(Cii+Cjj)/2-Cij` 的均值。C 默认为 Brownian 相关矩阵；`scale_vcv=false` 使用协方差，此时该差异等于 `patristic / 2`。`psr` 为 PSV 乘群落丰富度。两者直接使用树上统计量，方差也无需 VCV，接受稠密/CSC 输入和 `CommunityData`。

PSV 方差现按完整树物种池在指定丰富度下**均匀无放回**抽样定义，物种池包含未观测 tip；抽取完整物种池时方差为零。修复了与原 picante 相同的叶序依赖公式，原始协方差模式也计入不同 tip 方差的贡献。空群落和单物种群落返回 `NaN`，`compute_var=false` 时方差为 `NaN`；相关矩阵模式要求根到每个 tip 的长度为正。

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
diversity = psd(moment_data, tree)
diversity.pse         # [1/2, 5/6]：等丰度时 PSE = PSV
diversity.psc         # [1/2, 1/2]
pse(sparse([2.0 1 0 0]), tree).pse  # [4/9]
```

`pse` 按相对丰度加权物种对差异；`psc` 为各存在物种到其它物种的最小 `1-Cij` 的均值。原始协方差模式下 PSC 是一减最大共有枝长的均值，可为负，且依赖枝长单位；相关矩阵模式在非超度量树上选中的邻居可能与 patristic 最近物种不同。不同物种间的零长度路径保留。`psd` 一次返回四项指标、PSV/PSR 方差和丰富度，沿用单项调用的定义及 NaN 规则，接受并保留稠密/CSC 与 `CommunityData` 输入；算法及分配开销见[独立性能说明](performance_zh.md#解析多样性矩)。

稀疏数据的 PSV/PSR、PSE/PSC 使用观测物种连接树，保留原根对相关性和原始协方差的贡献。独立 PSV 只计算变异性，PSR 由同一路径得到；PSD 为每个群落只建一次连接树，供 PSV/PSE/PSC 共用。PSV/PSR（含 PSD 返回的对应项）的方差仍基于完整原树物种池，包括未观测物种及仅在单物种群落出现的物种。稠密场景保留顺序遍历，复杂度及实测见 [PSV/PSR 性能](performance_zh.md#psvpsr-连接树)和 [PSD 共用计算](performance_zh.md#psepsc-与-psd-共用计算)。

`expected_pd` 及别名 `pd_variance` 保留另一种模型：每个 tip 以 `richness/tree.ntips` 的概率**独立出现**，richness 是期望物种数。`upper_bound=false` 给出精确 Bernoulli 方差，默认 true 给出 picante 的协方差上界；并非固定样本量稀释曲线。枝长须有限、非负且总和有限；[算法与扩展性](performance_zh.md#解析多样性矩)单独说明。

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

两种模式都会在各群落合计观测物种较少时自动使用更小的连接树。先在原树上校验并按原丰度合并，再通过累加原始枝长压缩路径；覆盖树中大部分物种时保留整树路径，无需增加参数。

默认按丰度加权。`abundance_weighted=false` 只将各群落内部的存在物种改为等权；合并群落的组成和群落间权重仍使用原丰度，因此 beta/Fst 可能为负。这是现有 Julia 选项，Picante 没有对应参数。

空群落的 Dkk 及 Dkl/H 整行整列为 `NaN`；即使群落权重为零，也使 alpha/beta/Fst 为 `NaN`，但 total 仍可有效。全空输入的 total 为 `NaN`；零输入行返回空数组和 alpha=0。总多样性为零时 Fst 未定义，beta 也为零时返回 `NaN`。稠密/CSC 矩阵和 `CommunityData` 均遵循这些规则，输入不被修改；[算法与存储开销](performance_zh.md#rao-多样性分解)单独说明。

### 群落间系统发育距离

`comdist` 对两个群落间的全部物种对距离求均值；对角线包含同种的零距离及不同物种间距离。默认每个存在物种等权，`abundance_weighted=true` 使用丰度权重。空群落对应的整行、整列均为 `NaN`。

`comdistnt` 先分别求两个方向的最近目标物种距离均值。默认两个方向等权；`abundance_weighted=true` 时，方向内部按来源物种丰度加权，再按两个群落的原总丰度合并方向均值。没有合格目标的来源物种被跳过，其余权重重新归一化；仅一个方向有效时使用该方向，都无效时返回 `NaN`。默认非空群落的对角线为零。

```julia
distance_data = CommunityData(sparse([2.0 1 0; 0 1 3; 0 0 0]);
    sample_labels=["west", "east", "empty"], species_labels=["A", "B", "C"])
cross_mean = comdist(distance_data, tree; abundance_weighted=true)
nearest = comdistnt(distance_data, tree; exclude_conspecifics=true)
cross_mean[3, :] # [NaN, NaN, NaN]
nearest[1, 1]   # 2.0：第一个群落内部到不同物种的最近距离
requested_pairs = [(1, 2), (2, 1), (1, 1), (3, 2)]
selected_mean = comdist(distance_data, tree;
    abundance_weighted=true, pairs=requested_pairs)
# [10/3, 10/3, 8/9, NaN]
selected_nearest = comdistnt(distance_data, tree;
    abundance_weighted=true, exclude_conspecifics=true, pairs=requested_pairs)
# [22/7, 22/7, 2, NaN]
```

`comdist` 和 `comdistnt` 均支持 `pairs`，按请求顺序返回距离向量，保留重复、反向和自身对；索引始终对应完整输入的群落行序，空请求返回空向量。即使某群落不在请求中，也会校验全部数值、标签和树枝长。计算仅保留请求端点，重复/反向请求复用一次求值；`comdistnt` 为每个唯一对计算双向均值，同一目标的遍历可复用。使用 `mean_dist` 汇总非对角群落对时，传入包含 `distances`、`sample_labels` 和 `pairs` 的命名元组，并显式指定输出 `pairs`；`mean_dist` 要求自身距离为零或 NaN，而 `comdist` 的自身对包含群落内多样性。`comdistnt` 默认自身对为零，开启 `exclude_conspecifics=true` 后可能非零。

`exclude_conspecifics=true` 按物种身份排除同种；零长枝上的不同物种仍可具有零距离，此时单物种群落的对角线为 `NaN`。这些规则沿用现有 Julia 行为。Picante 在丰富度或总丰度不同的情况下采用不同的方向合并方式，排除同种时会排除全部零距离，并返回对角线为零的 `dist` 对象。

两者接受稠密/CSC 矩阵或 `CommunityData`，保持群落顺序，要求枝长有限、非负且总和有限。输入不会修改；矩阵可通过 `tip_labels` 描述列序，`CommunityData` 则按其标签对齐。CSC 路径无需稠密群落或完整物种距离矩阵；默认返回完整群落矩阵，指定群落对时仅返回请求值。[算法与存储开销](performance_zh.md#群落间系统发育距离)单独说明。

重复计算两种距离时，可以准备共用快照：

```julia
distance_context = prepare_community_distances(distance_data, tree)
reused_mean = comdist(distance_context; abundance_weighted=true,
    pairs=requested_pairs)
reused_nearest = comdistnt(distance_context; abundance_weighted=true,
    exclude_conspecifics=true, pairs=requested_pairs)

fixed_context = prepare_community_distances(distance_data, tree;
    pairs=requested_pairs)
fixed_nearest = comdistnt(fixed_context; abundance_weighted=true,
    exclude_conspecifics=true) # 请求顺序与 selected_nearest 相同
```

`CommunityDistanceContext` 保留丰度和存在两套权重、一份独立树拓扑及两种距离的工作区。准备时不传 `pairs`，后续可计算完整矩阵或更换群落对，共用覆盖全部已准备群落的拓扑。准备时固定 `pairs`，则只保留有效非空端点，并缓存请求调度；后续默认返回该请求顺序的向量，更换范围或顺序需要重新准备。

准备时校验全部输入数值、物种标签、枝长及群落标签，群落标签须唯一且非空。不保留原矩阵或树，之后修改输入不会改变快照；分析变化后的数据需重新准备。上下文字段属于内部实现，不应修改；返回数组各自独立，同一上下文的调用通过锁依次使用工作区，不缓存距离矩阵。用于 PD/PE/beta 的存在型 `CommunityContext` 保留另一组汇总，继续使用原入口。

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
