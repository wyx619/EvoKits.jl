# EvoKits 文档

[项目首页](../README_zh.md) · [English](README.md)

按分析任务选择指南。每页提供独立示例输入，页内代码块按顺序运行。

## 从这里开始

群落分析建议先阅读[输入与标签](community_input_zh.md)，再选择[多样性或 beta 指标](community_zh.md)，最后用[分区与排序](regionalization_zh.md)解释完整距离矩阵。[项目快速示例](../README_zh.md#快速开始)展示了这条分析路径。

| 分析任务 | 指南 |
| --- | --- |
| 读写、查询、修剪或补全树，读取/抽样后验树，时间切片 | [树操作](trees_zh.md) |
| 导入矩阵/长表，对齐物种，导出数据，选择复用上下文 | [群落输入](community_input_zh.md) |
| PD、特有性、多样性矩、Rao 分解，两两/多群落 beta | [群落多样性](community_zh.md) |
| 平均/最近分类单元距离、分支重叠、指定群落对 | [群落距离](community_distances_zh.md) |
| 选择零模型，计算 SES、beta SES、NRI/NTI | [零模型与群落结构](null_models_zh.md) |
| 性状汇总、物种聚合、功能分组和谱系分组 | [性状与聚合](aggregation_zh.md) |
| 聚类、固定 k GoM、后验平均距离、NMDS 及坐标对齐 | [分区与排序](regionalization_zh.md) |
| 性状的 lambda、K/PIC、D 信号 | [系统发育信号](signals_zh.md) |
| EDGE/EDGE2、热点/冷点和指标数值分级 | [保护指标](conservation_zh.md) |

NRI/NTI 描述相对于物种池的群落结构；lambda、K、D 描述性状在系统发育树上的信号，两者分别说明。

## 查阅与开发

| 查阅内容 | 入口 |
| --- | --- |
| 排序、结果形状、字段及未定义值 | [API 约定](api_zh.md) |
| 计算路径、复杂度、内存和基准条件 | [性能说明](performance_zh.md) |
| 实现文件及其加载依赖 | [源码组织](development_zh.md) |

指南解释模型并展示用法，API 汇总调用约定，性能文档集中说明算法、复杂度和实测。函数 docstring 可用 Julia 帮助模式查询，例如 `?expected_pd`。
