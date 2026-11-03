# 保护指标

[首页](../README_zh.md) · [English](conservation.md)

EDGE/EDGE2 得分、灭绝概率、指标热点/冷点筛选与数值分级。

以下示例使用这棵四物种树：

```julia
using EvoKits

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
```

## EDGE

`edge` 将完整树的 fair-proportion ED 与类别对应的 GE 得分结合：`log1p(ED) + GE*log(2)`。默认映射遵循 phyloregion 实际使用的概率值：

| 类别 | LC | NT | VU | EN | CR | EW/EX | DD |
|---|---|---|---|---|---|---|---|
| GE | 0.001 | 0.01 | 0.1 | 0.67 | 0.999 | 1 | missing |

```julia
edge_result = edge(tree, ["LC", "NT", "VU", "CR"])
edge_result.species_labels
edge_result.ed
edge_result.edge
edge_result.status

# 需要整数类别编码约定时，显式提供 GE 映射。
integer_edge = edge(tree, ["LC", "NT", "VU", "CR"];
    ge_mapping = Dict("LC" => 0, "NT" => 1, "VU" => 2, "EN" => 3, "CR" => 4))
```

概率映射与常见整数类别编码不同，比较结果时须统一定义。自定义字典替换整个默认映射，须覆盖全部传入类别。`EDGEResult` 按完整树 tip 顺序返回标签、ED、GE、EDGE、status 和实际映射。向量默认按树顺序，也可用 `tip_labels` 指定精确标签；`DataFrame` 默认读取 `species`/`category` 列，允许只覆盖部分物种。未提供记录、类别缺失或 DD 默认返回 missing 得分，状态分别为 `:absent_record`、`:missing_category`、`:unassessed`；`missing_category=:error` 改为报错。未知物种、未知类别与重复标签始终报错；枝长须有限非负。下方 EDGE2 使用独立模型和结果。

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

## 热点与冷点

`hotspots` 选择大于等于上尾分位数的值，`coldspots` 选择小于等于下尾分位数的值；`prob` 为尾部百分比，默认 2.5。两者采用 R type-7 分位数插值，保留阈值上的全部并列及原输入顺序和标签，实际选中比例可能超过请求比例。

```julia
metric = [0.0, 1.0, 2.0, 10.0, 12.0]
high = hotspots(metric; prob = 20, sample_labels = ["a", "b", "c", "d", "e"])
low = coldspots(metric; prob = 20)
high.threshold
high.selected

# 独立参考向量可让多份数据使用相同阈值。
shared_threshold = hotspots([1.0, 8.0, NaN]; reference = metric,
    prob = 20, invalid = :skip)
shared_threshold.valid_count
shared_threshold.reference_count
```

`hotspots(pd_result)`、`hotspots(pe_result)`、`hotspots(we_result)` 可直接接相应指标结果并保留群落标签，`coldspots` 同样支持。返回的 `MetricSelectionResult` 包括原值、筛选标记、阈值、尾部、百分比及有效次数。`invalid=:skip` 在 missing/NaN 原位置保留 missing 标记，Inf 始终报错；阈值参考中至少需一个有效值。零和负数正常参与。筛选直接使用分位数不等式，修正 R `coldspots` 在零阈值下因临时修改向量而误选高值的行为。

## 指标数值分级

`classify_values` 提供精确 Jenks 自然断点（默认 `:jenks`）、等距（`:equal`）和 type-7 分位数（`:quantile`）分级。

```julia
classes = classify_values(metric; k = 3, method = :jenks)
classes.classes       # 原输入顺序，类别从 1 开始
classes.upper_bounds
classes.breaks
classes.objective     # Jenks 类内平方和
classes.actual_k
```

`ValueClassificationResult` 还保留标签、原值、请求类别数、方法和原因。类别包含其上界，相同值始终同类；`breaks` 为最小值加全部上界，最小值单独成类时前两个断点相同。Jenks 类别数不超过不同值数量；等距/分位数方法先以有效观测数限制区间数，再去除重复上界。`actual_k` 统计上界数，等距/分位数区间可能为空；少于请求 k 时为 `:reduced_classes`，常量仅一类。空输入或全部无效时报错；`invalid=:skip` 在 missing/NaN 原位置保留 missing 类别，Inf 报错。

Jenks 在 m 个不同排序值上使用保留重复次数权重的精确动态规划，最小化类内平方和，耗时 O(k m²)、工作空间 O(k m)。默认 `max_work=50_000_000` 限制精确计算量，超出时报错。平移缩放和 Welford 区间方差避免原始平方和相消；原单位 SSE 可能溢出或下溢，`normalized_objective` 和 `objective_scale` 保留缩放结果。同代价时递归优先较晚的最后切分，可能与 R 断点不同而最优 SSE 相同。零值、负数、小数均保留；地图几何、颜色和绘图独立处理。
