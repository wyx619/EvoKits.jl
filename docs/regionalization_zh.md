# 生物地理分区与排序

[首页](../README_zh.md) · [English](regionalization.md)

对 beta 距离或群落–物种网络分区、评估结果、汇总后验树距离，进行 NMDS 排序并筛选区域指示物种。

## 基于 beta 多样性的生物地理分区

`phyloregion` 接收完整距离矩阵，默认使用 UPGMA（平均连接）层次聚类，按显式指定的 `k` 划分区域。也可直接接收完整的 `BetaDissResult`、`PhylobetaResult`，保留群落标签；`component` 可选 `:total`（默认）、`:turnover`、`:nestedness`。

```julia
using EvoKits, SparseArrays

regional_data = CommunityData(sparse([1 1 0 0; 1 1 1 0; 0 0 1 1; 0 1 1 1]);
    species_labels = ["A", "B", "C", "D"],
    sample_labels = ["west_1", "west_2", "east_1", "east_2"])
beta = beta_diss(regional_data)
regions = phyloregion(beta; k = 2)
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

排序坐标的计算见本页 [NMDS](#非度量多维尺度分析nmds)。

## 基于群落–物种网络的 Infomap

`infomap` 直接从群落矩阵或 `CommunityData` 构建无向网络：群落和物种分别作为两类节点，每个正值连接一条边。算法通过多次搜索最小化映射方程目标，产生分区数，无需完整 beta 距离矩阵。

```julia
network_regions = infomap(regional_data; seed = 42, trials = 10)
network_regions.membership       # 原始群落行序的分区归属
network_regions.k
network_regions.cluster_sizes
network_regions.codelength       # 完整群落–物种分区的编码长度，单位 bit

plain_network = infomap(regional_data; seed = 42, regularized = false)
weighted_network = infomap(regional_data; seed = 42, weights = :abundance)
```

默认开启正则化，`regularization_strength=1.0`。它改变流量模型，可能改变归属和分区数，甚至合并互不连通的部分。解释小型或稀疏网络时，可比较不同设置；不能凭不同流量模型下编码长度的大小判断分区优劣。这是在二部网络上运行普通 Infomap，其先验不施加专门的二部模型约束，也不保证地理连续。

本地 igraph 1.0.0 对照中，默认设置将两块各 2×2 的断开网络和 Africa 数据都合成一区。Africa 关闭正则化为三区，强度 0.1 在三个测试种子下均为四区。这些是敏感性结果，0.1 并非通用推荐强度。

`weights=:binary` 将任意正值视为存在，`:abundance` 使用正的有限 Float64 丰度作为边权。输入必须有限且非负。空群落报错；未出现物种和显式存储的零值不产生节点或边。CSC 输入不转换成稠密群落或邻接矩阵。谱系或功能分区可传入 `collapse_range`、`collapse_traits` 返回的 `.community`。

`InfomapResult` 按首个输入群落为区域编号，`cluster_sizes` 只统计群落，并记录边权模式、尝试次数、种子、正则化设置及底层 igraph 版本。`codelength` 对应包含物种节点的完整网络分区。相同输入顺序、设置和引擎版本下可通过种子复现，且不改变 Julia 的随机数状态；种子、输入重排或引擎版本变化可能产生不同局部最优解。

phyloregion 参考使用的旧 igraph 核心不支持正则化，其流量定义也与当前核心不同。因此关闭正则化可用于比较，但不保证分区或编码长度与旧 R 逐值一致。

## 区域指示物种

`indicators` 通过区域之间的有向 Poisson 或 Bernoulli KL 散度给物种评分。输入为外部模型提供的“物种 × 区域”theta 矩阵：`:poisson` 使用有限非负强度，`:bernoulli` 使用 [0,1] 概率。按原值计算，不自动归一化列。

```julia
theta = [0.9 0.1 0.2; 0.2 0.8 0.3; 0.1 0.3 0.7; 0.4 0.4 0.4]
indicator_result = indicators(theta; species_labels = ["A", "B", "C", "D"],
    region_labels = ["west", "east", "south"], method = :bernoulli,
    comparison = :min, top_indicators = 3, shared = false)
indicator_result.top_features  # 区域 × 请求排名，候选不足补 missing
indicator_result.scores
indicator_result.counts        # 每个区域实际返回的物种数
```

目标区域参数 p 与其他区域参数 q 的 Poisson KL 为 `p*log(p/q)+q-p`，Bernoulli KL 再加上 `1-p` 对 `1-q` 的对应项。`comparison=:min`（默认）取与其他区域的最小散度；`:max` 取最大散度。输入至少需要两个区域。

`shared=false` 时，物种只可入选其 theta 首个最大值所属的区域；`true` 允许同一物种出现在多个区域。分数并列按输入物种顺序排列；候选不足补 `missing`，不会重复填充物种。返回的 `IndicatorSpeciesResult` 保留标签、物种索引与名称、分数及参数。

零和一使用连续极限处理，散度可以为 Inf，零对零为零。计算逐区域复用一个物种评分向量，不保存“物种 × 区域 × 区域”临时张量。此接口不拟合 grade-of-membership 模型，也不从硬分区 membership 自动生成 theta；应先明确统计模型，再传入其估计值。

## 多棵后验树的平均距离

`mean_dist` 按群落标签汇总多棵树的 beta 距离，返回 `MeanDistResult`。输入迭代器只遍历一次，可边计算每棵树的 beta 边汇总。所有 beta 输入须为同一种结果、同一 Sørensen/Jaccard 指数族；`component` 选择 `:total`（默认）、`:turnover` 或 `:nestedness`。

```julia
posterior_trees = mktempdir() do dir
    path = joinpath(dir, "posterior.nwk")
    map(["((A:1,B:1):1,(C:1,D:1):1);",
         "((A:1,B:1):2,(C:1,D:1):2);"]) do newick
        write_tree(path, newick)
        read_tree(path)
    end
end

average = mean_dist((phylobeta(regional_data, tree) for tree in posterior_trees))
average.distances             # 对齐后的平均距离矩阵
average.counts                # 每个群落对的有限观测次数
average.replicates            # 实际消费的后验树数量
regions = phyloregion(average; k = 2)
context = prepare_phyloregion(average)
select_linkage(average)
optimal_phyloregion(context)

# 两个阶段都只处理指定群落对，保留重复、反向和自配对。
requested = [(1, 2), (2, 3), (2, 1)]
selected = mean_dist((phylobeta(regional_data, tree; pairs = requested)
    for tree in posterior_trees); pairs = requested)
selected.distances            # 按请求顺序返回向量
selected.counts
```

输出标签默认为首份输入的顺序，可用 `sample_labels` 显式指定完整集合。后续带标签输入可重排或缺少部分群落；未知标签、重名或空名报错。裸矩阵使用统一的输出标签顺序，维度须覆盖完整集合；各自带标签的矩阵写为 `(distances=D, sample_labels=labels)`，向量另提供自己的 `pairs`。普通距离与 beta 结果不能混合输入。

默认 `nan_policy=:skip` 对每个群落对分别跳过 NaN，只累计有限值。没有有限值时，均值为 NaN、次数为零，空群落或零 PD 群落的未定义对角线也保留。`:propagate` 遇到任意已观测 NaN 就返回 NaN，`counts` 仍报告有限次数；缺少群落或群落对不算一次观测。Inf、负距离、非零有限对角线和不对称矩阵报错。不同 beta 成分的有效次数可能不同，因此分别求均值后不保证仍满足原分解恒等式。

输出 `pairs` 索引输出标签，各份输入的 pairs 索引各自标签。输入中的重复或反向对须数值一致，每棵树只贡献一次；输出保留重复顺序。向量输入必须显式指定输出 pairs。指定对路径按标签和对数使用线性存储；完整输入矩阵仍会被完整校验，完整输出需要二次空间。输入不修改，缩放与补偿累加避免直接求和溢出并保留小量级距离。

完整、有限、对角线为零的平均结果可直接接入分区、聚类树复用、连接法比较和 k 评估；含未定义群落时需显式选择子集。该流程汇总树之间的距离，不计算区域归属的后验支持率。已完成 42 组本地 R 对照，其中 18 组来自真实多树 beta；普通测试只读取一张含四个小案例的固定 TSV，无需 R。

## 非度量多维尺度分析（NMDS）

`nmds` 将完整距离矩阵转换为排序坐标，通过单调回归和归一化 SMACOF 迭代优化全局 NMDS。目标为 Kruskal stress-1：`sqrt(sum((distance-disparity)^2)/sum(distance^2))`，按 0–1 尺度返回。也可直接接收完整 beta 和 `MeanDistResult`。

```julia
using Random

ordination = nmds(beta; component = :total, dims = 2, nstart = 20,
    rng = MersenneTwister(42))
ordination.points             # 群落 × 维度，保留输入行顺序。
ordination.stress
ordination.converged
ordination.reason
ordination.start_stress       # 每个起点的最终 stress。
ordination.stress_history     # 最佳起点的初始及各次接受后的 stress。

# PhyloregionResult 对区域排序；二维排序至少需要三个区域。
regional_ordination = nmds(phyloregion(beta; k = 3);
    nstart = 20, rng = MersenneTwister(42))
```

beta 输入选择 `:total`、`:turnover` 或 `:nestedness`，功能 beta 使用内部的 beta 结果。`nmds(average)` 使用 mean_dist 已选定的成分。`nmds(regions)` 使用区域间距离，标签 `"1"`、`"2"` 等对应 membership 的区域 ID；`nmds(context::PhyloregionContext)` 则对原始群落距离排序，不重建聚类树。

`nstart` 包含第一个起点，默认使用 PCoA，也可指定 `init=:random` 或“群落 × 维度”坐标矩阵；其余起点随机生成，由 `rng` 控制序列。选择最终 stress 最低的结果。每个起点默认 `maxiter=1000`、`tol=1e-7`（相对 stress 改善）和 `stress_tol=1e-4`；`maxiter=0` 仅评估起点。结果保留各起点的迭代数、收敛标记和终止原因，数值收敛和多起点均不保证全局最优。秩不足的起点可能始终停留在原子空间，随机起点可探索其它方向。

默认 `ties=:weak` 允许相等输入距离具有不同拟合值，与 vegan 默认并列含义一致；`:strong` 要求其拟合值相等。仅按精确相等分组，不使用 vegan 的绝对近似容差。零距离参与回归，不强制坐标重合。全为相等正值的距离在弱并列下可让任意坐标得到零 stress；全零输入则报错。要求 `1 <= dims < 群落数`，有解释力的 NMDS 通常需要远多于维度数的群落。NaN/Inf、非对称、非零对角线及指定对结果均报错，由调用者显式处理未定义群落。

`NMDSResult` 独立保存居中、到质心 RMS 距离为一且旋转至主轴的坐标。轴方向、相同主轴对应的方向及旋转不具有固定生态学含义，输入距离单位也不决定坐标单位。`dissimilarities`、`distances`、`disparities` 分别为原输入、坐标欧氏距离、单调拟合距离，按上三角逐列 `(1,2), (1,3), (2,3), (1,4), ...` 对齐，可用于 Shepard 图。不保留“群落对 × 迭代”历史。群落对存储为二次空间；默认 PCoA 初始化另需三次时间，可用 `init=:random` 避开该分解。当前入口不包含群落变换、物种坐标、step-across、half-change 缩放及地图/颜色生成。
