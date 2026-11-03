# EvoKits.jl

[English](README.md)

Julia 系统发育分析工具包，涵盖树操作、群落多样性、生物地理分区和比较分析，支持稠密与稀疏群落输入。

## 安装

```julia
using Pkg
Pkg.add(url = "https://github.com/wyx619/EvoKits.jl")
```

## 快速开始

创建一棵小树和带标签的群落数据，计算系统发育多样性（PD）、系统发育特有性（PE），并划分为两个区域：

```julia
using EvoKits, SparseArrays

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,C:2);")
    read_tree(path)
end

data = CommunityData(
    sparse([1 0 1; 0 1 1; 1 1 0; 1 1 1]);
    sample_labels = ["north", "south", "east", "west"],
    species_labels = ["A", "B", "C"],
)

pd(data, tree).pd
phylo_endemism(data, tree).pe

beta = phylobeta(data, tree)
regions = phyloregion(beta; k = 2)
regions.membership
```

矩阵以行为群落、列为物种；`CommunityData` 按标签将物种与树 tip 对齐。

## 使用指南

| 主题 | 内容 |
| --- | --- |
| [树操作](docs/trees_zh.md) | Newick 读写、拓扑、修剪、分类辅助补全、R/ape 映射及时间切片 |
| [群落多样性与结构](docs/community_zh.md) | 稀疏输入、PD、WE/PE、相对 PD、beta 分解、SES、NRI/NTI 及上下文复用 |
| [性状与聚合](docs/aggregation_zh.md) | 性状汇总、k-prototypes、组数选择、功能 beta 及谱系聚合 |
| [分区与排序](docs/regionalization_zh.md) | 八种连接法、混合归属、后验平均距离、NMDS 及 KL 指示物种 |
| [系统发育信号](docs/signals_zh.md) | Pagel's lambda、K、独立对比及 D |
| [保护指标](docs/conservation_zh.md) | EDGE/EDGE2、灭绝概率、热点/冷点及 Jenks 分级 |
| [API 约定](docs/api_zh.md) | 组合示例、排序、结果类型及边界情况 |
| [性能](docs/performance_zh.md) | 运行后端与性能基准 |
| [源码组织](docs/development_zh.md) | 内部目录职责与依赖顺序 |
