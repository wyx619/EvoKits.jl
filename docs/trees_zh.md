# 树操作

[首页](../README_zh.md) · [English](trees.md)

树的读写、查询、修剪、分类辅助补全、身份映射，以及不同时间深度的枝长切片。

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

加权距离使用枝长；设置 `weighted = false` 可改为计算边数。`distance`、`pairwise_distance` 和 `distances_to_tip` 在根深度相减会丢失精度时，累加原始分支段，保留长共有枝下的短距离。树源码按查询、操作、映射和指标分组，文件位置及依赖顺序见[源码组织](development_zh.md)。

## 后验树流式读取与抽样

`read_trees` 逐棵读取 Newick 或 NEXUS 文件，按各文件内容识别格式。Newick 记录以分号结束，可以跨行或多棵共用一行；支持方括号注释、根注释和单引号标签。可传入单个路径或按所需顺序排列的路径向量/元组。`burnin` 丢弃开头的树，`thin` 在剩余树中每隔指定棵数保留一棵，分别对每个文件应用。

```julia
using Random

posterior_example = mktempdir() do dir
    path = joinpath(dir, "posterior.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    write_tree(path, "((C:1,D:1):2,(B:1,A:1):2);"; append = true)
    write_tree(path, "((A:1,C:1):1,(B:1,D:1):1);"; append = true)
    summaries = read_trees(path; burnin = 1) do trees
        [evol_distinct(tree).scores for tree in trees]
    end
    sampled = tree_sampler(path; n = 2, replace = false,
        rng = MersenneTwister(42))
    (; summaries, sampled)
end
posterior_example.sampled.sample_indices
posterior_example.sampled.source_indices
```

返回的 `TreeStream` 只消费一次。读完或出错会关闭文件；do-block 在提前返回时也会关闭。手动管理时需在提前停止后调用 `close`。空文件不产生树。保留的树必须具有非空且唯一的 tip 标签；解析报错包含路径、原始树序号和行号。

`tree_sampler` 返回可迭代的 `TreeSampleResult`，默认 `n=100`、有放回，与 R 抽样器一致；`replace=false` 为无放回。每棵候选树等概率，文件树数不同时也不改变这一点。显式传入带种子的 `rng` 可复现 Julia 抽样，但不承诺与 R 的随机序列相同。候选总体为空、无放回数量超过总体均报错。

`trees`、`sample_indices`、`source_paths`、`source_indices` 和 `source_lines` 均按抽取顺序排列。抽样索引指向跨文件筛选后的总体，来源序号指向各文件的原始树号。重复抽取的树拥有独立数组。抽样扫描两遍，仅解析抽中的树；未抽中的记录仍检查分号、引号和注释是否闭合，不检查树结构。读取或抽样期间须保持文件内容不变。

### NEXUS 后验树文件

开头的 `#NEXUS` 标记启用 NEXUS 解析，不依赖扩展名。按文件顺序读取所有 `TREES` 块中的 `TREE`、`UTREE` 定义；`TRANSLATE` 只替换其所在块的叶标签，数字内部节点标签及根标签会保留。树名称及方括号注释会丢弃，包括 `[&R]`、`[&U]`；不根据这些标记重新定根。

```julia
nexus_example = mktempdir() do dir
    path = joinpath(dir, "posterior.trees")
    write(path, """
    #NEXUS
    BEGIN TREES;
      TRANSLATE 1 A, 2 B, 3 C;
      TREE STATE_0 = [&R] ((1:1,2:2):3,3:4);
      TREE STATE_100 = [&R] ((3:3,1:2):1,2:5);
    END;
    """)
    data = CommunityData(sparse([1 1 0; 1 0 1; 0 1 1]);
        species_labels = ["C", "A", "B"])
    average = read_trees(path) do trees
        mean_dist(data, trees)
    end
    sampled = tree_sampler(path; n = 2, replace = false,
        rng = MersenneTwister(42))
    (; average, sampled, regions = phyloregion(average; k = 2))
end
```

筛选与 `source_indices` 对同一文件跨块累计树号，不计元数据命令；`source_lines` 指向 `TREE`/`UTREE` 命令行。可混合传入 Newick 和 NEXUS 路径。每块至多一份翻译表，必须位于树定义前，键与物种名均须唯一；保留树的每个叶标签都须匹配表键。没有表时按字面标签读取，不自动改名，也不解析 `TAXA` 编号或 `LINK`（后者明确报错）。其他块跳过，每块必须用 `END` 或 `ENDBLOCK` 结束。burn-in 和抽样计数阶段仍校验命令与翻译表；只对实际解析的树校验 Newick 结构和叶标签映射。

两种迭代器均可接入[后验平均距离流程](regionalization_zh.md#多棵后验树的平均距离)。

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

## 分类信息辅助的树补全

`phylobuilder` 显式使用原树属级锚点或命名分类节点挂接缺失物种，返回 `PhyloBuilderResult`，包括独立的 `tree` 和每个请求物种的来源记录。

```julia
using DataFrames

source = mktempdir() do dir
    path = joinpath(dir, "source.nwk")
    write_tree(path, "((G_a:1,G_b:1)G:1,(H_a:1,I_a:1)Family:1)Root;")
    read_tree(path)
end
taxonomy = DataFrame(species = ["G_a", "G_c", "J_a", "Z_a"],
    genus = ["G", "G", "J", "Z"], family = [missing, missing, "Family", missing])
built = phylobuilder(taxonomy, source; extract = false)
built.status             # [:tree, :genus, :family, :missing]
built.added_species      # ["G_c", "J_a"]
built.unmatched_species  # ["Z_a"]
built.tree
```

表格的 species 列使用完整且精确的 tip 标签，不是种加词。也可输入标签向量，只有标签包含下划线或空格时才推断属。表格可选 genus 及更高分类列；默认 `taxonomy_cols=(:family,:suborder,:order,:class,:phylum)` 决定向上查找顺序，不自动修改标签拼写。

缺失物种挂接在原树同属 tip 的 MRCA；只有一个同属 tip 时挂其父节点。没有同属 tip 时先查同名内部属节点，再查更高分类节点。所有锚点取自原树，新增物种不能成为后续匹配依据。属非单系默认报错，需显式 `nonmonophyly=:mrca` 才允许挂其 MRCA；匹配的内部节点重名也报错。未能放置的物种在结果中保留为 `:missing`。

默认新 tip 枝长为锚点最长后代路径，原树须等时（`atol=0`、`rtol=1e-8`）。其他树须显式设置 `branch_length=:zero` 或有限非负末端枝长。原有边长保持；默认 `extract=true` 随后修剪到成功匹配的请求物种，可能收缩单子路径并移除节点标签；`false` 保留全部原 tip。输出重建节点编号和缓存，`anchor_nodes` 始终记录原树的紧凑节点 ID。`align_community` 继续拒绝未知物种；需要补全时先明确执行此步骤。

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
