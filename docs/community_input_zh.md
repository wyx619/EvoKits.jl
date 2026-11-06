# 群落输入与对齐

[文档目录](README_zh.md) · [首页](../README_zh.md) · [English](community_input.md)

`CommunityData` 将数值和标签保存在一起。矩阵以行为群落、列为物种，值须为有限非负实数；正值表示存在，丰度加权函数使用其大小。

## 矩阵、长表和索引三元组

```julia
using EvoKits, SparseArrays, DataFrames

data = CommunityData(sparse([2.0 1 0; 0 3 1; 0 0 0]);
    sample_labels = ["west", "east", "empty"], species_labels = ["A", "B", "C"])

records = DataFrame(sample = ["west", "west", "east", "empty"],
    species = ["A", "C", "A", "A"], value = [2, 1, 3, 0])
from_table = CommunityData(records; sample = :sample, species = :species, value = :value)
from_indices = CommunityData([1, 1, 2], [1, 2, 1], [2, 1, 3];
    sample_labels = ["west", "east"], species_labels = ["A", "C"])
```

| 输入 | 存储与标签 |
| --- | --- |
| 矩阵 | 复制数值，保留稠密/CSC 存储及数值类型；`storage=:sparse` 请求 CSC。 |
| 长 `DataFrame` | 直接构建 CSC，默认按标签首次出现顺序排列；`value=nothing` 按记录计数。 |
| 索引三元组 | 直接构建 CSC；提供完整 `sample_labels` 和 `species_labels`。 |

显式标签列表可保留空群落、全零物种列及指定顺序。标签须为唯一非空字符串或整数标识符，整数转为字符串。矩阵的默认标签为 `sample_1`、… 和 `species_1`、…；使用树时应提供实际物种名。缺失或非有限数值报错。

长表和三元组默认对重复记录求和；`duplicates=:error` 拒绝重复，包含零值记录，求和溢出也报错。整数丰度保留为计数，小数权重保留原值，Bool 记录转为 Int。

## 导出长表与重建群落

`community_longtable` 每个正值单元格导出一行，先按物种列序、再按群落行序排列。须保留其两份完整标签列表：正值记录无法单独恢复空群落和全零列。丰度不会展开为重复出现记录。

```julia
exported = community_longtable(data)
exported.table                     # sample、species、value 三列
restored = CommunityData(exported.table;
    sample_labels = exported.sample_labels, species_labels = exported.species_labels)
restored.matrix == data.matrix      # true
restored.sample_labels == data.sample_labels   # true
restored.species_labels == data.species_labels # true
```

零维和全空表同样适用。导出数组与输入独立；重建采用 CSC，并沿用长表数值规则。[导出开销](performance_zh.md#群落长表导出)单独说明。

## 树对齐与指标调用

```julia
tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
aligned = align_community(data, tree)
aligned.species_labels             # ["A", "B", "C", "D"]
pd(data, tree).pd                   # 自动按标签对齐
phylo_endemism(data, tree).pe
```

`CommunityData` 将缺失树物种补为零列，保留完整树和根。未知输入物种始终报错，即使其列全零；`align_community` 的 `missing_tips=:error` 还要求每个树物种均有列。对齐保留稠密/CSC 存储。在构造时设置标签，后续指标调用不要再传 `tip_labels`/`sample_labels`。

直接矩阵调用要求每个树 tip 均有一列，包含未观测物种。省略 `tip_labels` 时按 `tree.tip_labels` 解释列序；提供时，该关键字描述输入列序。结果保留输入群落行序。

## 重复分析

| 重复工作 | 预处理入口 | 复用内容 |
| --- | --- | --- |
| 仅匹配物种 | `align_community(data, tree)` | 对齐后的数值与标签 |
| PD、WE/PE、RPD、UniFrac/PhyloSor、beta | `prepare_community(data, tree)` | 分支集合、占据统计、群落汇总 |
| `comdist` / `comdistnt` 及两种权重 | `prepare_community_distances(data, tree)` | 丰度/存在记录、工作拓扑和距离工作区 |
| 同一聚类树的不同切分 | `prepare_phyloregion(beta)` | 距离快照和聚类树 |

```julia
ctx = prepare_community(data, tree)
pd(ctx).pd
phylo_endemism(ctx).pe
phylobeta(ctx; pairs = [(1, 2)])
```

构造、对齐和上下文拥有独立数组。数据及上下文字段应视为只读；原输入变化不会更新快照，须重新准备。指标结果可独立修改。各上下文保存的汇总不同，应选择函数支持的入口；它们不会缓存所有结果或完整两两矩阵。

后续可阅读[多样性与 beta](community_zh.md)、[群落距离](community_distances_zh.md)或[零模型](null_models_zh.md)。[API 约定](api_zh.md)汇总结果形状。
