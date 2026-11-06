# 生物地理分区与排序

[文档目录](README_zh.md) · [首页](../README_zh.md) · [English](regionalization.md)

本页：[聚类](#基于-beta-多样性的生物地理分区) · [评估分区](#复用聚类树及评估分区) · [GoM](#区域混合归属) · [指示物种](#区域指示物种) · [后验均值](#多棵后验树的平均距离) · [NMDS](#非度量多维尺度分析nmds) · [对齐与配色](#坐标对齐与配色)

对 beta 距离分区、拟合混合归属、评估结果、汇总后验树距离，进行 NMDS 排序并筛选区域指示物种。

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

固定行序下的并列规则可重复：新链从最早活动簇开始，最近距离相同时优先前驱，否则取最早活动簇；同高度合并保留发现顺序。质心/中位连接选活动槽编号字典序最早的一对。并列选择可能与 R 不同，也可能随行重排变化。输入不修改，结果数组独立；成本见[性能说明](performance_zh.md#分区与排序)。

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

`PhyloregionContext` 保存原始距离的压缩快照、标签和完整聚类树，修改原输入不影响快照；上下文内部数组应只读。重复切分时只需准备一次；`cutree` 仅返回归属，`phyloregion(ctx; k)` 还计算区域汇总。

`select_linkage` 也可接完整距离矩阵或 beta 结果，可用 `methods=(:average,:complete,:single)` 限定方法。它比较原距离与各群落对首次合并高度的相关性，分数越大表示距离保留越好，不代表生态解释必然正确。常量距离/高度和仅一个群落对的相关性为 NaN；全部未定义时 `best_method=nothing`。不构建完整 cophenetic 矩阵；上下文输入对相同方法复用已有树。

`optimal_phyloregion` 也可接完整距离矩阵或 beta 结果，只构建一次树；`k` 是评估上限，默认 `min(n,20)`。簇内 WSS 为该簇无序群落对的距离平方和除以群落数，TSS 对全体群落计算，BSS=TSS−WSS，`ev=BSS/TSS`。对非欧氏差异，这些是距离评分，EV 可能不单调或超出 [0,1]。

建议 k 遵循本地 phyloregion 的阈值规则：先按降序尝试 `ev_thresholds=(0.95,0.9,0.8,0.75,0.67,0.5,0.33)`，再按升序尝试 `inc_thresholds=(0.01,0.05,0.1)`。分别找首个达到解释率阈值的 k、首个下一步增量足够小的 k，对首个成功阈值组合推荐两者较大值；容差为 `10^(-precision)`，默认 precision=3，评估末端 k 无下一步增量。这是阈值启发式。阈值不满足或全零离差时返回 `optimal_k=nothing` 及完整表，`reason` 分别为 `:thresholds_not_met`、`:zero_dispersion`。平方和按原始距离平方单位输出，可能溢出/下溢；EV 和判断使用归一化距离避免该溢出。

排序坐标的计算见本页 [NMDS](#非度量多维尺度分析nmds)。

## 区域混合归属

`fitgom` 接收群落矩阵或 `CommunityData`，拟合固定 `k` 个潜在区域。一个群落可以同时具有多个区域的组成特征：`omega` 给出混合权重，`theta` 给出各区域的物种组成概率。权重表达组成混合程度，不是分类正确的置信度。

```julia
using StableRNGs

gom = fitgom(regional_data; k = 2, nstart = 5, rng = StableRNG(42))
gom.omega                 # 群落 × 区域，每行和为一
gom.theta                 # 物种 × 区域，每列和为一
gom.membership            # 每个群落首个最大权重所属区域
gom.converged, gom.reason, gom.residual
gom.start_objectives      # 所有起点，包括未收敛的起点
gom_indicators = indicators(gom; top_indicators = 3)
```

输入按原值使用：二值矩阵拟合存在型组成，整数丰度对应多项式模型，分数权重对应加权似然。标量 `shape` 默认为 `1/(k×物种数)`。最大化的 NEF 目标为 `sum(x*log(q)) + shape*sum(log(theta)) + sum(log(omega))/k`，其中 `q[i,j] = sum(omega[i,r]*theta[j,r] for r in 1:k)`。`shape` 直接作为正伪计数使用，不减一；整体缩放计数会改变先验与数据的权衡。

优化采用 EM，只接受使目标进一步上升的对数比外推；`accelerate=false` 使用纯 EM，每轮均含两步 EM。`init_theta`、`init_omega` 可提供第一个起点的正初值，按列/行归一化；其余起点随机初始化，返回目标最大的解。`history`、`start_*`、`maxiter`、`tol`、`accelerate` 保留过程和设置。`residual` 是再做一步普通 EM 时两个因子的最大绝对变化，只有 `residual <= tol`（默认 `1e-8`）才标记收敛；终止原因为 `:tolerance`、`:maxiter` 或 `:stalled`。仍可能得到局部最优解。目标对齐 maptpx 的标量 shape 版本，其条件 SQP/EM 优化路径与这里不同。

群落和物种保留输入顺序，区域按 `omega` 总权重降序排列，精确并列保留原序；硬归属的 `cluster_sizes` 可以包含零。空群落报错，全零物种列保留，参与先验的物种数并在 `theta` 中保留正的先验贡献。输入不修改；存储见[性能说明](performance_zh.md#分区与排序)。谱系/功能模型传入聚合后的 `.community`。`1 <= k <= min(群落数,物种数)`；自动选择 k、Bayes factor 和不确定性估计保留在长远规划。

## 区域指示物种

`indicators` 通过区域之间的有向 Poisson 或 Bernoulli KL 散度给物种评分。可直接传入 `GoMResult`，或外部模型提供的“物种 × 区域”theta 矩阵：`:poisson` 使用有限非负强度，`:bernoulli` 使用 [0,1] 概率。按原值计算，不自动归一化列。

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

零和一使用连续极限处理，散度可以为 Inf，零对零为零。计算逐区域复用一个物种评分向量，不保存“物种 × 区域 × 区域”临时张量。评分器不从硬分区 membership 自动估计 theta；可先运行 `fitgom` 或传入其他统计模型的估计值。

## 多棵后验树的平均距离

`mean_dist(data, trees)` 逐棵计算并直接累计一个系统发育 beta 成分，返回 `MeanDistResult`。群落记录与请求调度只准备一次，每棵树分别匹配和校验，不再为每棵树创建三组成分结果。树迭代器只消费一次；`index_family` 选择 `:sorensen` 或 `:jaccard`，`component` 选择 `:total`（默认）、`:turnover` 或 `:nestedness`。原有 `mean_dist(inputs)` 继续接受已经计算的 beta 或距离结果迭代器。

```julia
using Random

requested = [(1, 2), (2, 3), (2, 1)]
posterior_analysis = mktempdir() do dir
    path = joinpath(dir, "posterior.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    write_tree(path, "((C:1,D:1):2,(B:1,A:1):2);"; append = true)
    average = read_trees(path) do trees
        mean_dist(regional_data, trees)
    end
    # 只准备请求端点，每个唯一群落对计算一次。
    selected = read_trees(path) do trees
        mean_dist(regional_data, trees; pairs = requested)
    end
    sampled = tree_sampler(path; n = 4, rng = MersenneTwister(42))
    sampled_average = mean_dist(regional_data, sampled)
    (; average, selected, sampled, sampled_average)
end

average = posterior_analysis.average
average.distances             # 对齐后的平均距离矩阵
average.counts                # 每个群落对的有限观测次数
average.replicates            # 实际消费的后验树数量
regions = phyloregion(average; k = 2)
context = prepare_phyloregion(average)
select_linkage(average)
optimal_phyloregion(context)

selected = posterior_analysis.selected
selected.distances            # 按请求顺序返回向量
selected.counts
```

逐文件 burn-in、间隔保留和抽样来源见[树流与抽样](trees_zh.md#后验树流式读取与抽样)。`CommunityData` 按物种标签分别对齐每棵树，后验树的 tip 顺序可以不同。

相同调用也接受 NEXUS 后验树文件，并应用其中的 `TRANSLATE` 翻译表。

直接累计时，`CommunityData` 提供固定的群落与物种标签，可以省略树上的未观测物种。矩阵列须覆盖每棵树的完整 tip 集合；`tip_labels` 固定各列含义，省略时由第一棵树的 tip 顺序定义后续所有树的列含义。所有输入物种（包括零列）都须存在于每棵树；未请求的群落也完整校验数值。`block_size=256` 控制临时列块，完整输出仍需群落数的平方空间，详见[性能说明](performance_zh.md#均值累计存储)。

对于 `mean_dist(inputs)`，输出标签默认为首份输入的顺序，可用 `sample_labels` 显式指定完整集合。后续带标签输入可重排或缺少部分群落；未知标签、重名或空名报错。裸矩阵使用统一的输出标签顺序，维度须覆盖完整集合；各自带标签的矩阵写为 `(distances=D, sample_labels=labels)`，向量另提供自己的 `pairs`。普通距离与 beta 结果不能混合输入；所有 beta 结果须为同一种类和同一 Sørensen/Jaccard 指数族。

默认 `nan_policy=:skip` 对每个群落对分别跳过 NaN，只累计有限值。没有有限值时，均值为 NaN、次数为零，空群落或零 PD 群落的未定义对角线也保留。`:propagate` 遇到任意已观测 NaN 就返回 NaN，`counts` 仍报告有限次数；缺少群落或群落对不算一次观测。Inf、负距离、非零有限对角线和不对称矩阵报错。不同 beta 成分的有效次数可能不同，因此分别求均值后不保证仍满足原分解恒等式。

输出 `pairs` 索引输出标签，各份输入的 pairs 索引各自标签。输入中的重复或反向对须数值一致，每棵树只贡献一次；输出保留重复顺序。向量输入必须显式指定输出 pairs。指定对路径按标签和对数使用线性存储；完整输入矩阵仍会被完整校验，完整输出需要二次空间。输入不修改，缩放与补偿累加避免直接求和溢出并保留小量级距离。

完整、有限、对角线为零的平均结果可直接接入分区、聚类树复用、连接法比较和 k 评估；含未定义群落时需显式选择子集。该流程汇总树之间的距离，不计算区域归属的后验支持率。

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

`NMDSResult` 独立保存居中、到质心 RMS 距离为一且旋转至主轴的坐标。轴方向、相同主轴对应的方向及旋转不具有固定生态学含义，输入距离单位也不决定坐标单位。`dissimilarities`、`distances`、`disparities` 分别为原输入、坐标欧氏距离、单调拟合距离，按上三角逐列 `(1,2), (1,3), (2,3), (1,4), ...` 对齐，可用于 Shepard 图。存储与初始化成本见[性能说明](performance_zh.md#分区与排序)。当前入口不包含群落变换、物种坐标、step-across 及 half-change 缩放。

## 坐标对齐与配色

`align_nmds` 按共同标签，通过正交 Procrustes 最小二乘将已有坐标对齐到参考。默认只平移、旋转和镜像，保持源坐标距离；`scale=true` 另估计一个统一尺度，`reflection=false` 限制为不含镜像的旋转。独立返回坐标与诊断，保留原 NMDS 拟合及 stress。

```julia
reference_colors = nmds_colors(ordination)
comparison = nmds(beta; dims = 2, nstart = 5, rng = MersenneTwister(43))
aligned = align_nmds(comparison, ordination)
aligned.matched_labels
aligned.rmse                      # 参考坐标单位下的锚点距离均方根。

comparison_colors = nmds_colors(aligned; scale = reference_colors.scale,
    outside = :clip)
comparison_colors.hex             # #RRGGBB，保留 comparison 的原始行序。
comparison_colors.rgb             # [0,1] 范围的量化 sRGB 通道。
comparison_colors.clipped         # 超出固定参考范围的行。

regional_colors = nmds_colors(regional_ordination)
site_colors = regional_colors.hex[phyloregion(beta; k = 3).membership]
```

对齐也可接坐标矩阵，通过 `sample_labels`、`reference_labels` 提供标签，省略时使用 `sample_1` 等默认标签。支持行重排与部分共有群落；共同群落决定变换，再应用到源坐标的全部行，保留原序。两份输入须有相同的 `d` 个轴，至少 `d+1` 个共同标签，居中锚点及其交叉乘积须满列秩；不足或秩退化时报错。`source_indices`、`reference_indices`、`rotation`、`scale`、`source_center`、`reference_center` 保存变换，`residuals`、`rmse` 描述共同锚点。不同分区结果的区域编号不自动表示同一生态实体；跨分析应使用稳定群落标签或先明确区域对应关系。

`nmds_colors` 接 NMDS/对齐结果或单轴、双轴坐标矩阵，可提供矩阵的 `sample_labels`。固定色板对齐 phyloregion 的 `hexcols`：HCL 边界色、中间 sRGB 舍入、每轴 200 级索引。两轴统一使用参考中的最大跨度，保留长宽比例；颜色距离不保证等于生态差异。恒定轴映射到色板中点，单轴输入补零；单点或全部恒定的坐标获得统一中心色。只计算请求的颜色。

多次比较时先对齐至共同参考，再复用其 `.scale`，保持范围和单位。仅复用范围不会对齐方向。默认拒绝超界坐标，显式 `outside=:clip` 才截断颜色坐标，并保留受影响行的标记；实际坐标不变。色板网格边界可能受浮点舍入影响。输入不修改，无 GIS、绘图或新增 Julia 依赖，输出颜色可交给绘图工具使用。
