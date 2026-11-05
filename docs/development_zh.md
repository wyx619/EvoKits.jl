# 源码组织

[首页](../README_zh.md) · [English](development.md)

EvoKits 保持单一 Julia 模块。[EvoKits.jl](../src/EvoKits.jl) 声明公共导出，并按依赖顺序加载源码。子目录用于按职责组织同一模块内的文件，公共名称与调用方式保持一致。

入口文件按功能分组导出，结果类型紧邻对应函数。NRI/NTI 和 SES 归入空模型与群落结构；系统发育信号和保护指标各自独立成组。源码加载区按依赖阶段标注，一个阶段可以跨越多个目录。MKL 保持为首个加载的依赖。

| 目录 | 职责 | 代表源码 |
|---|---|---|
| `base/` | 共用标签、树表示、读写及数值辅助 | [labels.jl](../src/base/labels.jl) |
| `tree/queries/` | 校验、查询、距离与谱系分组 | [clade_groups.jl](../src/tree/queries/clade_groups.jl) |
| `tree/operations/` | 修剪、变换、多分叉处理、时间切片与补全 | [phylobuilder.jl](../src/tree/operations/phylobuilder.jl) |
| `tree/mapping/` | 节点/边身份映射及 R/ape 编号 | [mapping.jl](../src/tree/mapping/mapping.jl) |
| `tree/metrics/` | 进化独特性与树形指标 | [metrics.jl](../src/tree/metrics/metrics.jl) |
| `traits/` | 性状对齐与系统发育信号 | [alignment.jl](../src/traits/alignment.jl) |
| `community/core/` | 群落输入、长表导出、分支遍历、上下文与指标适配 | [input.jl](../src/community/core/input.jl)、[longtable.jl](../src/community/core/longtable.jl) |
| `community/diversity/` | PD、MPD/MNTD、群落间距离、特有性、Rao 及 beta 多样性 | [community_distances.jl](../src/community/diversity/community_distances.jl)、[rao.jl](../src/community/diversity/rao.jl)、[beta.jl](../src/community/diversity/beta.jl) |
| `community/aggregation/` | 物种/谱系聚合、功能分组与组数选择 | [functional_beta.jl](../src/community/aggregation/functional_beta.jl) |
| `community/regionalization/` | 分区、连接法评估、平均距离、NMDS 与指示物种 | [regionalization.jl](../src/community/regionalization/regionalization.jl) |
| `community/null_models/` | 随机化、SES、beta SES 与 NRI/NTI | [ses.jl](../src/community/null_models/ses.jl) |
| `community/conservation/` | EDGE/EDGE2、热点/冷点与数值分级 | [conservation.jl](../src/community/conservation/conservation.jl) |

共用标签校验位于 `base/labels.jl`：将非空字符串或整数标识符转换为唯一字符串标签，沿用已有默认标签与报错规则。树补全与群落输入使用同一套辅助函数。base 和 tree 先于 traits、community 加载，因此 `phylobuilder` 可在群落定义加载前独立使用。群落文件之间的依赖决定加载顺序，同一子目录的文件不要求连续加载。

`base/tree_stream.jl` 切分记录并解析 Newick 树体，`base/tree_nexus.jl` 处理 NEXUS 块、命令和叶标签翻译，`tree/operations/tree_sampling.jl` 通过两遍扫描抽样，仅解析抽中的记录。这些文件均先于 traits 和 community 加载，不依赖群落层。

`gom.jl` 提供固定 k 混合归属结果及指示物种适配，`gom_fit.jl` 实现稀疏似然/EM 累计和有目标上升保护的对数比加速，使用现有 Julia 依赖。

`nmds_alignment.jl` 按共同标签计算正交 Procrustes 变换；`nmds_colors.jl` 将坐标映射到 phyloregion 固定色板，并保存可复用范围。两者均使用现有依赖，在独立坐标副本上操作。

`rao.jl` 定义 Rao 熵及其分解，复用 `community_distances.jl` 准备的稀疏记录与分支距离工作区；两者均先于 `data_metrics.jl` 的 `CommunityData` 适配加载。

[comdist_virtual.jl](../src/community/diversity/comdist_virtual.jl) 为 `comdist` 选择整树/压缩树路径并调度指定群落对，一次后序遍历压缩未观测及单分叉路径、累加原始枝长，随后复用共用的两遍距离内核。[comdistnt_virtual.jl](../src/community/diversity/comdistnt_virtual.jl) 复用该连接树，按目标组织双向请求，重复/反向对去重后再合并方向均值。Rao 内部路径继续使用整树遍历。

[distance_context.jl](../src/community/core/distance_context.jl) 在 `CommunityDistanceContext` 中拥有丰度/存在记录、最小距离拓扑及两种工作区，固定对准备仅保留有效端点并缓存双向调度；锁保护共享工作区，结果数组始终独立。原 `CommunityContext` 继续保留存在型分支汇总，不增加距离存储。

`phylogenetic_species.jl` 用稀疏记录与树统计量计算 PSV/PSR/PSE/PSC/PSD，包括有限物种池方差及 PSD 共用准备，并复用 `community_distances.jl` 的不同目标身份工作区。`expected_pd.jl` 定义 Bernoulli 模型下的解析 PD 矩与分支丰度；两者先于保留 PD 变体及 bootstrap 的 `diversity.jl` 加载，随后由 `data_metrics.jl` 添加输入适配。

数值测试继续放在 `test/`，小型离线参考表放在 `test/fixtures/`。结构调整采用相关测试文件和包加载检查。本地验收脚本、历史源码快照保存在受 Git 忽略的 `validation/`；历史 `git show` 命令中的路径对应原提交的目录结构。
