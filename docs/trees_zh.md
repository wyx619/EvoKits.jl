# 树操作

[首页](../README_zh.md) · [English](trees.md)

树的读写、查询、修剪、身份映射，以及不同时间深度的枝长切片。

## 树的读写与拓扑

```julia
using EvoKits, SparseArrays

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end

ultrametric_tree = make_ultrametric(tree)
binary_tree = multi2di(ultrametric_tree; random = false)
write_tree("binary.nwk", binary_tree)
```

`CompactTree` 是 EvoKits 工具的通用输入类型。它缓存拓扑、枝长、遍历顺序和 tip 标签查找信息。`read_tree` 将 Newick 文件直接读取为这种表示。

树查询明确使用紧凑节点 ID、tip 位置或 tip 标签：

```julia
mrca(tree, "A", "B")
distance(tree, 1, 2)
pairwise_distance(tree)
descendant_tips(tree, tree.root; labels = true)
tree_imbalance(tree; metric = :sackin)
```

加权距离使用枝长；设置 `weighted = false` 可改为计算边数。源代码分为 `base`、`tree`、`traits` 和 `community` 四层。

## 树的修剪与映射

```julia
pruned = keep_tip(tree, ["A", "C", "D"])
without_b = drop_tip(tree, "B")

ref = build_phyloref(pruned)
node_map = phyloref_node_table(pruned; map = ref)
edge_map = phyloref_edge_table(pruned; map = ref)

ape_edges = ape_edge_matrix(pruned; map = ref)
signature = phylo_edge_signature(pruned, 1; map = ref)
```

`PhyloRef` 提供稳定的节点和边身份标识，包括 `CompactTree` ID 与 R/ape 编号之间的转换。映射表和签名方便在不同工具及运行环境之间比较树结果。

## 树枝长时间切片

`timeslice(tree, depth; collapse=false)` 从定年等时树中移除距今 `depth` 时间内的枝长，返回新的 `CompactTree`，保留全部物种及原列序。旧分支保持原长，跨截面的分支截短，起点恰好位于截面或更年轻的分支置零。depth 使用树枝长单位，合法范围为零至树高；零保持枝长，树高使全部分支归零。

```julia
mktempdir() do dir
    path = joinpath(dir, "tree.nwk")
    write_tree(path, "(((A:1,B:1):1,C:2):2,(D:2,E:2):2);")
    tree = read_tree(path)
    data = CommunityData(sparse([1 1 0 0 0; 0 0 1 0 0; 0 0 0 1 1]);
        species_labels = tree.tip_labels)
    profiles = map((0.0, 1.0, 2.0)) do depth
        sliced = timeslice(tree, depth; collapse = true)
        ctx = prepare_community(data, sliced)
        (; depth, pe = phylo_endemism(ctx).pe,
            regions = phyloregion(phylobeta(ctx); k = 2))
    end
end
```

默认保留拓扑、节点编号、边序及标签，并重建根距离缓存。`collapse=true` 折叠所有恰为零的内部边（包括原有零边），保留根及全部 tip，重新编号；被折叠节点的标签随节点移除。任何正枝长都保留。含根路径的 PD、PE 和按枝长计算的 beta 在折叠前后保持一致；边序诊断与拓扑统计可能变化。

验证沿用 `get_clades` 的 `atol=0.0`、`rtol=1e-8` 与最长后代路径年龄。容差内近似等时树按各父节点的局部年龄切片；容差不移动截面，也不决定折叠边。时间与额外空间为 O(节点数+物种数)，输入不修改。

比较多个绝对时间深度时，各次都从原树切片，并为每棵切片树重新准备群落上下文。切到树高时，所有 PD/PE 为零，beta 为 NaN，因此该端点不能直接用于分区。
