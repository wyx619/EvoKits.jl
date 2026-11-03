# 源码组织

[首页](../README_zh.md) · [English](development.md)

EvoKits 保持单一 Julia 模块。[EvoKits.jl](../src/EvoKits.jl) 声明公共导出，并按依赖顺序加载源码。子目录用于按职责组织同一模块内的文件，公共名称与调用方式保持一致。

| 目录 | 职责 | 代表源码 |
|---|---|---|
| `base/` | 共用标签、树表示、读写及数值辅助 | [labels.jl](../src/base/labels.jl) |
| `tree/queries/` | 校验、查询、距离与谱系分组 | [clade_groups.jl](../src/tree/queries/clade_groups.jl) |
| `tree/operations/` | 修剪、变换、多分叉处理、时间切片与补全 | [phylobuilder.jl](../src/tree/operations/phylobuilder.jl) |
| `tree/mapping/` | 节点/边身份映射及 R/ape 编号 | [mapping.jl](../src/tree/mapping/mapping.jl) |
| `tree/metrics/` | 进化独特性与树形指标 | [metrics.jl](../src/tree/metrics/metrics.jl) |
| `traits/` | 性状对齐与系统发育信号 | [alignment.jl](../src/traits/alignment.jl) |
| `community/core/` | 群落输入、分支遍历、上下文与指标适配 | [input.jl](../src/community/core/input.jl) |
| `community/diversity/` | PD、MPD/MNTD、特有性及 beta 多样性 | [beta.jl](../src/community/diversity/beta.jl) |
| `community/aggregation/` | 物种/谱系聚合、功能分组与组数选择 | [functional_beta.jl](../src/community/aggregation/functional_beta.jl) |
| `community/regionalization/` | 分区、Infomap、连接法评估、平均距离、NMDS 与指示物种 | [regionalization.jl](../src/community/regionalization/regionalization.jl) |
| `community/null_models/` | 随机化、SES、beta SES 与 NRI/NTI | [ses.jl](../src/community/null_models/ses.jl) |
| `community/conservation/` | EDGE/EDGE2、热点/冷点与数值分级 | [conservation.jl](../src/community/conservation/conservation.jl) |

共用标签校验位于 `base/labels.jl`：将非空字符串或整数标识符转换为唯一字符串标签，沿用已有默认标签与报错规则。树补全与群落输入使用同一套辅助函数。base 和 tree 先于 traits、community 加载，因此 `phylobuilder` 可在群落定义加载前独立使用。群落文件之间的依赖决定加载顺序，同一子目录的文件不要求连续加载。

`infomap.jl` 负责群落校验及群落分区结果，`infomap_backend.jl` 隔离 IGraphs 原生对象、随机状态及调用。EvoKits 直接依赖 IGraphs，Graphs 作为其间接依赖引入。

数值测试继续放在 `test/`，小型离线参考表放在 `test/fixtures/`。结构调整采用相关测试文件和包加载检查。本地验收脚本、历史源码快照保存在受 Git 忽略的 `validation/`；历史 `git show` 命令中的路径对应原提交的目录结构。
