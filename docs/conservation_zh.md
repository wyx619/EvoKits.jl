# 保护指标

[首页](../README_zh.md) · [English](conservation.md)

EDGE2 得分与物种灭绝概率工具。

以下示例使用这棵四物种树：

```julia
using EvoKits

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
```

## EDGE2

EDGE2 将进化独特性与物种灭绝概率结合，使用与其他工具相同的 `CompactTree`，无需额外的树解析器或 `Phylo.jl` 依赖。

```julia
pext = [0.10, 0.20, 0.30, 0.40]  # 按 tree.tip_labels 排序
edge2 = edge2_mod(tree, pext)

edge2.species
edge2.edge       # EDGE2 得分
edge2.ed         # 进化独特性
edge2.pd         # 系统发育多样性总量
edge2.epd_loss   # 预期 PD 损失
```

参考实现中的字段拼写 `EDGE`、`ED`、`TBL`、`PD` 和 `ePDloss` 仍可作为兼容属性使用。灭绝概率必须是 `(0, 1]` 内的有限值，并按树 tip 顺序排列。

基于类别的工作流可使用 `PextDistribution`、`calc_pext_distribution`、`load_pext_distribution_csv` 和 `sample_pext`。`reorder_pext_to_tree` 会先验证物种覆盖情况，再按树顺序返回概率。
