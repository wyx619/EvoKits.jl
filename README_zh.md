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

`ses_pd`、`ses_mpd` 和 `ses_mntd` 只保留一份随机化结果及可复用的物种/树工作区，不再构建稠密“群落 × 物种”随机化缓冲区。支持 `:taxa_labels`、`:richness`、`:frequency`、`:sample_pool`、`:phylogeny_pool`、`:trialswap` 和 `:independentswap`。存储重构保留原抽样规则及基于排列的随机数消耗：richness 保留行和，frequency 保留列和，交换模型同时保留两者。现有 `:sample_pool` 与 `:richness` 一样使用完整树物种池，`:phylogeny_pool` 还会置换物种标签。丰度加权 MPD/MNTD 仍随机化存在状态，并将原丰度固定在树 tip 位置；零权重的随机物种仍作为 MNTD 最近邻候选。

三种 SES 均按有限的随机值计算均值和样本标准差：没有有效随机值时，均值、标准差和 z 均为 `NaN`；只有一个有效随机值时，标准差及 z 为 `NaN`；零模型标准差为零时，z 为 `NaN`。秩 p 值保留原约定 `(1 + 小于观测值的次数 + 0.5 × 相等次数) / (runs + 1)`，非有限随机值不增加秩计数。PD 已修复重复抽样间分支访问标记复用导致的漏计；不包含根时直接累加 MRCA 以下分支，避免用相减去除很长的共有根路径。

SES 虚拟边在深度相减可能丢失精度时，也会回退到原始分支累加，以保留长共享路径下的极短距离。

如需 Caper 兼容的多分叉 ED 处理方式，可使用 `evol_distinct(tree; type = :fair_proportion, polytomy_correction = :isaac)`；其他支持的模式为 `:mooers` 和 `:none`。

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

## 群落性状与物种分组

`map_trait` 按群落汇总一个数值性状，`aggregate_species` 按既定分组聚合物种。两者均接受稠密/CSC 矩阵或 `CommunityData`；性状和分组可用按物种列顺序的向量、按物种名索引的字典或 `DataFrame` 提供。按名称精确匹配，表格物种禁止重复。矩阵入口可指定 `species_labels` 和 `sample_labels`。

```julia
using EvoKits, SparseArrays, DataFrames

data = CommunityData(sparse([2 1 0; 0 3 1; 0 0 0]);
    species_labels = ["A", "B", "C"], sample_labels = ["AB", "BC", "empty"])
traits = DataFrame(species = ["C", "A", "B"], trait = [-2.0, 2.0, 6.0])
summary = map_trait(data, traits; statistic = :mean, abundance_weighted = true)
summary.values # [10/3, 4, NaN]
summary.matched_richness # [2, 2, 0]

groups = Dict("A" => "group_1", "B" => "group_1", "C" => "group_2")
grouped = aggregate_species(data, groups) # 各群落每组的物种数
proportions = aggregate_species(data, groups;
    abundance_weighted = true, normalize = true)
beta_diss(grouped) # 分组后的 CommunityData 可直接用于分类指标
```

`map_trait` 返回 `CommunityTraitResult`，包含汇总 `values`、群落标签、原始 `richness`、实际使用的 `matched_richness` 和 `weight_sum`，以及选项信息。`statistic` 支持 `:sum`（默认）、`:mean`、`:sd`。默认每个存在物种计一次；`abundance_weighted=true` 使用原丰度。加权标准差将丰度解释为可靠性权重，校正分母为 `sum(w) - sum(w^2)/sum(w)`；单位权重时等于样本标准差。`corrected=false` 使用总体标准差分母 `sum(w)`。加权均值和标准差不受权重整体缩放影响，整数丰度不会展开成重复记录。数值汇总使用 Float64，溢出时报错；求和补偿相消误差，标准差先中心化再在线更新。

默认遇到存在物种缺少性状或分组时报错；可设置 `missing_traits=:skip` 或 `missing_groups=:skip` 跳过。未出现的物种可以没有映射，多余映射物种会忽略。数值性状用 Julia `missing` 表示缺失；匹配到的 NaN、Inf 或非数值性状会报错。空群落和全部被跳过的群落保留标签：总和为零，均值/标准差为 `NaN`，使用权重为零。仅一个匹配物种时，校正标准差为 `NaN`，总体标准差为零。

`aggregate_species` 返回新的 `CommunityData`，CSC 列对应物种组。默认统计每组存在的物种数；`abundance_weighted=true` 汇总丰度并保留数值类型（Bool 转 Int）。`normalize=true` 返回各群落保留分组的 Float64 比例，空行仍为零。组名使用字符串或整数，默认按输入物种列首次出现的组排序，包含已映射零列对应的组；`group_labels` 可固定完整顺序并保留额外空组。性状表列由 `species`/`trait` 选择，分组表列由 `species`/`group` 选择。输入不会被修改，CSC 输入不展开为稠密群落矩阵或长出现记录。分组由调用者提供，物种聚类与地理绘图分别处理。

## 混合性状与功能 beta

`functional_groups` 根据数值与类别性状，通过 k-prototypes 将物种分为 `k` 组；`collapse_traits` 将群落按该映射聚合；`functional_beta` 对功能组的存在/缺失计算 Sørensen 或 Jaccard beta 分解。

```julia
using EvoKits, DataFrames, SparseArrays, Random

data = CommunityData(sparse([1 1 0 0; 0 0 1 1; 1 0 1 0]);
    species_labels = ["A", "B", "C", "D"], sample_labels = ["west", "east", "mixed"])
traits = DataFrame(species = ["C", "A", "D", "B"],
    height = [10.0, 0.0, 12.0, 2.0], form = ["b", "a", "b", "a"])

groups = functional_groups(data, traits; k = 2, numeric_scale = :none,
    gamma = 10, nstart = 10, rng = MersenneTwister(42))
groups.membership             # [1, 1, 2, 2]，按 data 的物种列序
counts = collapse_traits(data, groups; statistic = :richness)
fb = functional_beta(data, groups)
fb.community.matrix           # 稀疏功能组群落
fb.beta.total                 # 群落两两功能 beta 总差异
regions = phyloregion(fb.beta; k = 2)

# 也可直接拟合并计算。
fb = functional_beta(data, traits; k = 2, rng = MersenneTwister(42))
```

聚类距离为数值性状的加权平方差，加上 `gamma` 乘以类别性状的加权不匹配数。默认 `numeric_scale=:sd` 使用样本标准差；也可选择 `:range` 或 `:none`。常量列的尺度设为一。`trait_weights=Dict(:height => 2.0, :form => 0.5)` 可单独设置正权重；默认 `gamma=:auto` 令全部物种对的平均数值贡献与类别贡献相等，任一贡献为零时取一。类别值支持字符串、Symbol、Bool，以及通过 `categorical_traits` 显式指定的数值编码；同一数值类别列必须使用相同具体类型。

表格拟合使用全部物种；数据对象拟合按群落的完整物种列序选择，包括全零列。物种按名称精确匹配；缺失映射或性状默认报错，`missing_traits=:skip` 显式跳过并在结果的 `excluded_species` 中报告。`trait_cols` 可指定使用的性状列。

每次随机起点选取不同物种；多次起点返回目标值最低的结果。空组由距离最远且所在组多于一个物种的成员填充，类别众数并列取字典序较小值。结果提供 `objective`、`iterations` 和 `converged`；这是局部优化，不能保证全局最优。需要复现时传入固定状态的 `rng`。结果中的数值中心采用 `(值 - numeric_offsets) / numeric_scales` 坐标，物种映射可跨群落复用。

`collapse_traits` 支持 `:presence`、`:richness`、`:abundance` 及行比例归一化，保留全部功能组和空群落。`functional_beta` 可传 `pairs` 只计算指定群落对；其 beta 始终基于功能组是否出现，`statistic=:abundance` 只改变返回的聚合群落数值。未实现 R 的 `quick_elbow` 自动选组数；本实现的尺度和初始化规则是上面明确列出的约定，不承诺与 R 随机分组逐项相同。

已用本地 clustMixType 0.5.2 对照 40 个案例：对齐尺度、权重与初始原型后，分组、目标值、原型、稀疏聚合及 beta 一致。自动权重默认值、并列分配和重复/空组处理按各自约定；普通 Julia 测试读取六个固定 R 案例，不需要 R 环境。

## 树分组与谱系群落聚合

`get_clades` 将完整的定年等时树按“距今时间”切分为互不重叠的谱系组，返回 `CladeGroupingResult`，包括树 tip 顺序的 `tip_labels`、`membership`，按首个 tip 排序的 `group_labels`、`group_sizes`、`group_nodes`（CompactTree 原始节点 ID），以及截面元数据。树不修剪、不修改。

```julia
mktempdir() do dir
    path = joinpath(dir, "tree.nwk")
    write_tree(path, "(((A:1,B:1):1,C:2):2,(D:2,E:2):2);")
    tree = read_tree(path)
    groups = get_clades(tree; cut = 2)
    groups.membership     # [1,1,1,2,2]
    groups.group_labels   # ["clade_1", "clade_2"]
    groups.group_nodes    # 原始树节点 ID，不是 tip 列位置

    data = CommunityData(sparse([2 4 0 0 0; 0 0 5 1 2; 0 0 0 0 0]);
        species_labels = tree.tip_labels, sample_labels = ["ab", "cde", "empty"])
    collapsed = collapse_range(data, groups) # 复用完整树分组
    collapsed.community.matrix              # [1 0; 1 1; 0 0]
    counts = collapse_range(data, tree; cut = 2, statistic = :richness)
    counts.community.matrix                 # [2 0; 1 2; 0 0]
    mass = collapse_range(data, groups; statistic = :abundance)
    mass.community.matrix                   # [6 0; 5 3; 0 0]
    beta_diss(collapsed.community)           # 可直接接分类 beta
end
```

`cut` 和目标 `k` 必须且只能指定一项。显式截面上的节点，默认 `boundary=:older` 保留整支为一组，`:younger` 则展开。截面等于树高时 older 为一组；cut=0 时 older 可能保留零龄谱系，younger 则全部为单 tip。边界精确比较当前表示的年龄，不加 epsilon。

目标 `k` 按完整分支事件展开，多分叉或同龄节点一起处理，实际分组数可能超过目标。例如上树 `get_clades(tree; k=3)` 返回 `requested_k=3`、`actual_k=4`、`cut=2`、`boundary=:younger`。用结果的 cut 与 boundary 再调用可重现同一分组。k 模式不另指定 boundary，k=1 保留根，合法范围 1:ntips；单子节点不增加组数。

节点年龄由最长后代路径计算，避免长根干下的短近期枝被相减丢失。等时性比较根到 tip 的最长/最短路径，默认 `atol=0.0, rtol=1e-8`；容差仅校验树，不合并近似同龄事件。接受的近似等时树也按最长后代路径取节点年龄。枝长须有限且非负、路径和有限；非等时树报错，由调用者先明确完成定年或变换。显式截面耗时/存储 O(节点数+tip数)，目标 k 另排序分支事件，不构造完整 tip 距离矩阵或为每个节点存全套后代集合。

`collapse_range` 接树及截面参数，或接可复用的 `CladeGroupingResult`。返回 `LineageCommunityResult`，包含 CSC `CommunityData` 的 `community`、独立分组 `grouping`、`statistic` 与 `normalize`。默认 `:presence` 每群落每谱系至多计一次；`:richness` 计组内存在物种数；`:abundance` 汇总原丰度并保留类型（Bool 转 Int），检查求和溢出。`normalize=true` 在选定统计量之后归一化，空行仍为零。矩阵入口用 `tip_labels`、`sample_labels`，默认按完整树/分组的 tip 列序。

物种按名称精确匹配，未出现树物种贡献零值，完整树的全部组列保留；任何未知物种均报错，包括零列。不因物种缺失而稠密补零，也不转成长表或构造稠密分组矩阵。输入不修改，结果数组独立；可复用分组的内部数组应只读。聚合列已是谱系组，可接 `beta_diss`；原物种树不能直接当作这些组的树使用。本操作也不同于截短枝长的树时间切片。

## 树枝长时间切片

`timeslice(tree, depth; collapse=false)` 从定年等时树中移除距今 `depth` 时间内的枝长，返回新的 `CompactTree`，保留全部物种及原列序。旧分支保持原长，跨截面的分支截短，起点恰好位于截面或更年轻的分支置零。depth 使用树枝长单位，合法范围为零至树高；零保持枝长，树高使全部分支归零。

```julia
mktempdir() do dir
    path = joinpath(dir, "tree.nwk")
    write_tree(path, "(((A:1,B:1):1,C:2):2,(D:2,E:2):2);")
    tree = read_tree(path)
    data = CommunityData(sparse([1 1 0 0 0; 0 0 1 0 0; 0 0 0 1 1]);
        species_labels = tree.tip_labels)
    profiles = map((0.0, 1.0, 2.0)) do depth
        sliced = timeslice(tree, depth; collapse = true)
        ctx = prepare_community(data, sliced)
        (; depth, pe = phylo_endemism(ctx).pe,
            regions = phyloregion(phylobeta(ctx); k = 2))
    end
end
```

默认保留拓扑、节点编号、边序及标签，并重建根距离缓存。`collapse=true` 折叠所有恰为零的内部边（包括原有零边），保留根及全部 tip，重新编号；被折叠节点的标签随节点移除。任何正枝长都保留。含根路径的 PD、PE 和按枝长计算的 beta 在折叠前后保持一致；边序诊断与拓扑统计可能变化。

验证沿用 `get_clades` 的 `atol=0.0`、`rtol=1e-8` 与最长后代路径年龄。容差内近似等时树按各父节点的局部年龄切片；容差不移动截面，也不决定折叠边。时间与额外空间为 O(节点数+物种数)，输入不修改。

比较多个绝对时间深度时，各次都从原树切片，并为每棵切片树重新准备群落上下文。切到树高时，所有 PD/PE 为零，beta 为 NaN，因此该端点不能直接用于分区。

## 基于 beta 多样性的生物地理分区

`phyloregion` 接收完整距离矩阵，默认使用 UPGMA（平均连接）层次聚类，按显式指定的 `k` 划分区域。也可直接接收完整的 `BetaDissResult`、`PhylobetaResult`，保留群落标签；`component` 可选 `:total`（默认）、`:turnover`、`:nestedness`。

```julia
regional_data = CommunityData(sparse([1 1 0 0; 1 1 1 0; 0 0 1 1; 0 1 1 1]);
    species_labels = ["A", "B", "C", "D"],
    sample_labels = ["west_1", "west_2", "east_1", "east_2"])
regions = phyloregion(beta_diss(regional_data); k = 2)
regions.membership       # 原始群落行序，区域编号从 1 开始
regions.cluster_sizes
regions.region_distances
regions.distinctiveness

# 也可直接传完整距离矩阵及标签。
regions = phyloregion([0 1 8 10; 1 0 6 12; 8 6 0 2; 10 12 2 0];
    k = 2, sample_labels = ["A", "B", "C", "D"])
# membership = [1, 1, 2, 2]；region_distances = [0 9; 9 0]
```

结果 `PhyloregionResult` 包含群落标签及归属 `membership`、分区大小 `cluster_sizes`、区域间差异 `region_distances`、区域独特性 `distinctiveness`，以及完整聚类树 `merge`、`heights` 和选项 `k`、`method`、来源 `component`。区域按其首个输入群落的顺序编号。`merge` 负数表示原始输入行号，正数表示先前合并行号。区域间差异为两区域全部跨群落对的平均距离，对角线设零；区域独特性为它与其他区域平均差异的均值，每个区域权重相同，不按群落数加权。`k=1` 时独特性为 `NaN`。

距离输入必须非空、方形、严格对称、有限、非负且对角线为零，不要求欧氏距离。转换为 Float64 后须保持有限正距离。任何未定义值都报错，包括空群落造成的 NaN 对角线；需要调用者先显式选择定义有效的子集及其标签，不会自动删除群落或填补距离。指定群落对的 beta 结果不能替代完整距离输入。

支持以下 8 种方法，递推语义对齐 R `hclust`：

| `method` | 含义 |
|---|---|
| `:average` | UPGMA，全部跨群落距离的平均值 |
| `:single` | 单连接，最小跨群落距离 |
| `:complete` | 完全连接，最大跨群落距离 |
| `:mcquitty` | WPGMA，两个子簇等权更新 |
| `:ward_D` | R 旧版 Ward，直接递推输入距离 |
| `:ward_D2` | 平方距离上的 Ward 递推，输出平方根高度 |
| `:centroid` | UPGMC，在输入数值上应用质心递推 |
| `:median` | WPGMC，在输入数值上应用中位连接递推 |

`Symbol("ward.D")`、`Symbol("ward.D2")` 是接受的别名，结果统一记录为 `:ward_D`、`:ward_D2`。Ward 的最小方差解释要求欧氏距离。质心/中位连接不自动将输入平方，允许合并高度下降，保留实际合并顺序。更新距离超出有限 Float64 范围时报错，需先缩放输入。

六种单调方法使用 O(n²) 最近邻链，质心/中位连接使用全局最小值加行最近邻缓存，最坏时间 O(n³)。工作区均用 `n*(n-1)/2` 个 Float64 保存压缩三角；完整输入和区域输出仍需二次方内存。固定行序下的并列规则可重复：新链从最早活动簇开始，最近距离相同时优先前驱，否则取最早活动簇；同高度合并保留发现顺序。质心/中位连接选活动槽编号字典序最早的一对。并列选择可能与 R 不同，也可能随行重排变化。输入不修改，结果数组独立。

### 复用聚类树及评估分区

```julia
ctx = prepare_phyloregion(beta_diss(regional_data); method = :average)
cutree(ctx; k = 2)                  # 仅归属，O(n)，不重新聚类
regions3 = phyloregion(ctx; k = 3)  # 复用树，另算区域统计
cutree(regions3; k = 2)             # PhyloregionResult 也保留完整树

selection = select_linkage(ctx)     # 默认比较 8 种方法
selection.methods
selection.correlations             # 原距离与树上合并距离的 Pearson 相关
selection.best_method              # 首个有限最高分方法，或 nothing

evaluation = optimal_phyloregion(ctx; k = 4) # 评估 1 到 4 的全部切分
evaluation.df                      # k、ev、totwss、totbss、tss
evaluation.optimal_k                # 建议分区数，或 nothing
evaluation.reason
```

`PhyloregionContext` 保存原始距离的压缩快照、标签和完整聚类树，修改原输入不影响快照；上下文内部数组应只读。准备时保留一份压缩三角，另用一份临时三角聚类。重复 `cutree` 只切树；`phyloregion(ctx; k)` 仍会遍历群落对计算区域均值，并分配 k×k 区域矩阵。

`select_linkage` 也可接完整距离矩阵或 beta 结果，可用 `methods=(:average,:complete,:single)` 限定方法。它比较原距离与各群落对首次合并高度的相关性，分数越大表示距离保留越好，不代表生态解释必然正确。常量距离/高度和仅一个群落对的相关性为 NaN；全部未定义时 `best_method=nothing`。不构建完整 cophenetic 矩阵；上下文输入对相同方法复用已有树。

`optimal_phyloregion` 也可接完整距离矩阵或 beta 结果，只构建一次树；`k` 是评估上限，默认 `min(n,20)`。簇内 WSS 为该簇无序群落对的距离平方和除以群落数，TSS 对全体群落计算，BSS=TSS−WSS，`ev=BSS/TSS`。对非欧氏差异，这些是距离评分，EV 可能不单调或超出 [0,1]。

建议 k 遵循本地 phyloregion 的阈值规则：先按降序尝试 `ev_thresholds=(0.95,0.9,0.8,0.75,0.67,0.5,0.33)`，再按升序尝试 `inc_thresholds=(0.01,0.05,0.1)`。分别找首个达到解释率阈值的 k、首个下一步增量足够小的 k，对首个成功阈值组合推荐两者较大值；容差为 `10^(-precision)`，默认 precision=3，评估末端 k 无下一步增量。这是阈值启发式。阈值不满足或全零离差时返回 `optimal_k=nothing` 及完整表，`reason` 分别为 `:thresholds_not_met`、`:zero_dispersion`。平方和按原始距离平方单位输出，可能溢出/下溢；EV 和判断使用归一化距离避免该溢出。全部评分只遍历每个群落对一次，额外存储 O(n)。

本轮不新增依赖；NMDS、地理约束、空间多边形和绘图尚未实现。

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
- `weighted_endemism` 返回 `WeightedEndemismResult`，包含 `sample_labels`、`we`、`richness` 和按输入物种列顺序排列的 `species_occupancy`。
- `phylo_endemism` 返回 `PhyloEndemismResult`，包含 `sample_labels`、`pe`、`richness`、按树的边顺序排列的 `branch_occupancy` 和 `weighted`。
- `relative_pd` 返回 `RelativePDResult`，包含 `sample_labels`、`rpd`、`pd_normalized`、`pd_equal_normalized`、`richness`、`total_branch_length` 和 `positive_branch_count`。
- `phylobeta` 返回 `PhylobetaResult`，包含 `sample_labels`、`total`、`turnover`、`nestedness`、`index_family` 和 `pairs`。
- `ses_pd`、`ses_mpd` 和 `ses_mntd` 返回带类型的结果，包含观测值、零模型均值和标准差、z 分数、秩 p 值、运行次数和零模型。空群落的 MPD/MNTD 为 `NaN`；单物种群落也为 `NaN`，但正质量的丰度加权 MPD 为零。
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
