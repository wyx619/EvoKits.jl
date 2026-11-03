# 性状与群落聚合

[首页](../README_zh.md) · [English](aggregation.md)

汇总群落性状，按预定义分组、功能组或谱系组聚合物种。

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

`collapse_traits` 支持 `:presence`、`:richness`、`:abundance` 及行比例归一化，保留全部功能组和空群落。`functional_beta` 可传 `pairs` 只计算指定群落对；其 beta 始终基于功能组是否出现，`statistic=:abundance` 只改变返回的聚合群落数值。本实现的尺度和初始化规则是上面明确列出的约定，不承诺与 R 随机分组逐项相同。

已用本地 clustMixType 0.5.2 对照 40 个案例：对齐尺度、权重与初始原型后，分组、目标值、原型、稀疏聚合及 beta 一致。自动权重默认值、并列分配和重复/空组处理按各自约定；普通 Julia 测试读取六个固定 R 案例，不需要 R 环境。

### 功能组数评估与选择

`optimal_functional_groups` 拟合一组候选 k，并使用 phyloregion 的 `quick_elbow` 所依赖的 `lbow` 启发式给出建议。物种匹配、数值缩放、性状权重和 gamma 只准备一次，保持各候选目标值可比。

```julia
selection = optimal_functional_groups(data, traits; k = 2:4,
    nstart = 20, rng = MersenneTwister(42))
selection.k
selection.objectives
selection.converged
selection.reason
if selection.grouping !== nothing
    selected_beta = functional_beta(data, selection.grouping)
end
```

候选必须包含至少两个不同的正整数，且不超过性状完整的物种数。程序按 k 递增拟合，在候选和起点间顺次消耗传入的 RNG。其他拟合参数沿用 `functional_groups`，也支持只传性状表。返回的 `FunctionalGroupSelectionResult` 包括全部候选的目标值、迭代数、收敛标记，以及 `optimal_k` 和可直接复用的 `grouping`；结果只保留所选分组。

启发式在归一化目标曲线上使用 `low=0.08`、`max_fraction=0.9`。它不保证生态意义上的最优组数。全零目标、候选未收敛、目标曲线增加或阈值无法建议时，`reason` 分别为 `:zero_objective`、`:nonconverged_candidates`、`:nonmonotone_objective`、`:thresholds_not_met`，建议 k 和 grouping 均为 `nothing`，完整诊断仍保留。使用建议前应检查目标曲线和收敛情况。

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
