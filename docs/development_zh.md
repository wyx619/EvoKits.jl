# 源码组织

[文档目录](README_zh.md) · [首页](../README_zh.md) · [English](development.md)

EvoKits 保持单一 Julia 模块。[EvoKits.jl](../src/EvoKits.jl) 按功能分组导出公共名称，并按依赖顺序加载源码。子目录用于组织同一模块内的职责，不引入子模块。

## 查找实现

| `src/` 下的目录 | 职责 | 入口 |
| --- | --- | --- |
| `base/` | 共用标签、树表示、读写及数值辅助 | [labels.jl](../src/base/labels.jl)、[tree_io.jl](../src/base/tree_io.jl) |
| `tree/queries/` | 校验、查询、距离与谱系分组 | [clade_groups.jl](../src/tree/queries/clade_groups.jl) |
| `tree/operations/` | 修剪、变换、时间切片、补全与抽样 | [phylobuilder.jl](../src/tree/operations/phylobuilder.jl) |
| `tree/mapping/` | 节点/边身份及 R/ape 编号 | [mapping.jl](../src/tree/mapping/mapping.jl) |
| `tree/metrics/` | 进化独特性与树形指标 | [metrics.jl](../src/tree/metrics/metrics.jl) |
| `traits/` | 性状对齐与系统发育信号 | [alignment.jl](../src/traits/alignment.jl) |
| `community/core/` | 输入、长表导出、分支汇总、上下文与适配 | [input.jl](../src/community/core/input.jl)、[data_metrics.jl](../src/community/core/data_metrics.jl) |
| `community/diversity/` | PD、物种多样性、距离、特有性、Rao 及 beta | [community_distances.jl](../src/community/diversity/community_distances.jl)、[phylogenetic_species.jl](../src/community/diversity/phylogenetic_species.jl) |
| `community/aggregation/` | 性状/物种/谱系分组与功能 beta | [functional_beta.jl](../src/community/aggregation/functional_beta.jl) |
| `community/regionalization/` | 聚类、GoM、后验均值、NMDS 与指示物种 | [regionalization.jl](../src/community/regionalization/regionalization.jl) |
| `community/null_models/` | 随机化、SES、beta SES 与 NRI/NTI | [ses.jl](../src/community/null_models/ses.jl) |
| `community/conservation/` | EDGE/EDGE2、热点/冷点与数值分级 | [conservation.jl](../src/community/conservation/conservation.jl) |

## 加载与共用机制

MKL 先于其他依赖加载，base 和 tree 先于 traits、community。共用标签校验位于 `base/labels.jl`，树补全因此可独立于群落输入使用。加载阶段可以跨目录，同目录文件不要求连续加载。

| 跨功能共用内容 | 文件 |
| --- | --- |
| 树流切分、NEXUS 命令及选中树抽样 | [tree_stream.jl](../src/base/tree_stream.jl)、[tree_nexus.jl](../src/base/tree_nexus.jl)、[tree_sampling.jl](../src/tree/operations/tree_sampling.jl) |
| PD/PE/beta 的存在型分支集合与占据统计 | [branches.jl](../src/community/core/branches.jl)、[context.jl](../src/community/core/context.jl) |
| 稀疏记录、自动连接树与距离复用 | [community_distances.jl](../src/community/diversity/community_distances.jl)、[comdist_virtual.jl](../src/community/diversity/comdist_virtual.jl)、[distance_context.jl](../src/community/core/distance_context.jl) |
| PSV/PSE/PSC 及 PD 变体的群落树 | [species_virtual.jl](../src/community/diversity/species_virtual.jl)、[pd_workspace.jl](../src/community/diversity/pd_workspace.jl) |
| 既有距离与后验 beta 的均值累计 | [mean_dist.jl](../src/community/regionalization/mean_dist.jl)、[mean_beta.jl](../src/community/regionalization/mean_beta.jl) |

群落适配在调用的内核之后加载。Rao 复用距离记录/拓扑；PSD 共用群落准备；直接后验 beta 复用均值累计器，但逐树缓冲区独立于单树上下文。算法与复杂度见[性能说明](performance_zh.md)，快照/排序约定见 [API 约定](api_zh.md)。

## 修改与验证

数值测试放在 `test/`，小型离线 R 参考表放在 `test/fixtures/`，来源由各 fixture README 说明；普通测试不应要求安装 R 或联网。

行为改动运行相关模块测试。目录调整增加包加载及依赖模块检查；文档改动检查相对链接、中英文一致性和修改后的示例，源码未变时无需重跑完整数值测试。用法更新专题指南，约定更新 `api_zh.md`，成本更新 `performance_zh.md`。

本地验收脚本、历史快照和详细验收记录放在受忽略的 `validation/`，不作为公开 API 或文档依赖。历史 Git 路径对应原提交的目录结构。
