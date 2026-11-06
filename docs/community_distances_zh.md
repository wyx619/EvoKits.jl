# 群落间系统发育距离

[文档目录](README_zh.md) · [首页](../README_zh.md) · [English](community_distances.md)

按平均物种距离、最近分类单元距离或分支重叠比较群落。输入构造和标签对齐见[群落输入](community_input_zh.md)。

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

## 群落间系统发育距离

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


## UniFrac 与 PhyloSor

`unifrac(data, tree)` 返回非加权分支距离；`phylosor(data, tree)` 返回相似度 `2*shared_PD/(PD1+PD2)`。两者按存在信息计算、保留行序，返回完整矩阵。

UniFrac 在分支并集总长为零时返回零，零 PD 群落与正 PD 群落之间为一。PhyloSor 对正 PD 群落的对角线返回一，分母为零时返回零，包含空群落和零 PD 群落。零分母规则见[矩阵语义](api_zh.md#矩阵语义)。需要周转/嵌套分解或指定分支重叠群落对时，使用 `phylobeta`。

```julia
unifrac(data, tree)
phylosor(data, tree)
phylobeta(data, tree; pairs = [(1, 2)])
```

聚类要求完整、有限且对角线为零的距离矩阵。`comdist` 对角线包含群落内多样性，PhyloSor 则是相似度；作为距离使用前须明确解释和转换。[分区指南](regionalization_zh.md)说明了可接受的输入。
