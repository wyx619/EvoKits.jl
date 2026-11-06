# 零模型与群落结构

[文档目录](README_zh.md) · [首页](../README_zh.md) · [English](null_models.md)

本页：[群落内 SES](#群落内-ses) · [零模型选择](#选择零模型) · [NRI/NTI](#nri-和-nti) · [Beta SES](#系统发育-beta-ses)

SES 和 NRI/NTI 将群落结构与选定零模型比较，描述物种组合；性状的系统发育信号另见[系统发育信号](signals_zh.md)。

```julia
using EvoKits, SparseArrays

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
community = [1 0 1 0; 0 1 1 0]
data = CommunityData(sparse(community);
    sample_labels = ["west", "east"], species_labels = tree.tip_labels)
```

## 群落内 SES

```julia
using Random

pd_effect = ses_pd(data, tree; null_model = :richness,
    runs = 99, rng = MersenneTwister(42))
mpd_effect = ses_mpd(data, tree; null_model = :richness,
    runs = 99, rng = MersenneTwister(42))
mntd_effect = ses_mntd(data, tree; null_model = :richness,
    runs = 99, rng = MersenneTwister(42))
pd_effect.pd_obs_z
mpd_effect.mpd_nri             # mpd_obs_z 的相反数
mntd_effect.mntd_nti           # mntd_obs_z 的相反数
```

各结果保留输入群落行序、丰富度（`ntaxa`）、观测值、零模型均值和样本 SD、秩、z/p 值及运行设置；字段前缀分别为 `pd_`、`mpd_`、`mntd_`。表格输出使用 `ses_pd_dataframe`、`ses_mpd_dataframe` 或 `ses_mntd_dataframe`。

`z = (观测值 - 零模型均值)/零模型SD`。空群落 MPD/MNTD 为 `NaN`；单物种群落也为 `NaN`，但有正质量的加权 MPD 为零。PD 遵循[根与单物种规则](api_zh.md#边界情况)。丰度加权 MPD/MNTD 随机化存在状态，将原丰度固定在树 tip 位置；零权重随机物种仍为 MNTD 候选。

均值/SD 仅使用有限随机值：无有效值时均值/SD/z 为 `NaN`，只有一个时 SD/z 为 `NaN`，零 SD 也给出 `NaN` z。下尾中秩分数为 `(1 + 小于观测值的次数 + 0.5*相等次数)/(runs + 1)`；非有限值不增加秩计数。这不是双尾 p 值。

## 选择零模型

| 零模型 | 随机化及约束 |
| --- | --- |
| `:taxa_labels`（`:tipshuffle`） | 一次共用的树 tip 标签置换，保留分类学重叠关系。 |
| `:richness`、`:sample_pool` | 各群落独立从完整树物种池抽取原丰富度的物种。 |
| `:phylogeny_pool` | 丰富度抽样后再进行共用标签置换。 |
| `:frequency` | 保持物种占据群落数，群落丰富度可改变。 |
| `:trialswap`、`:independentswap` | 棋盘交换保持行列边际；每次重复从观测群落重新开始。 |

完整树物种池包含未观测物种；`:sample_pool` 与 `:richness` 使用同一物种池。交换模型由 `iterations` 和 `max_attempts_per_swap` 控制，有限交换次数不保证均匀抽样。比较同一抽样下的统计量时，使用相同种子和设置的新 RNG。零模型结果在线汇总，不保留全部随机群落。

## NRI 和 NTI

NRI 和 NTI 描述相对于零模型的群落系统发育结构，分别是 MPD 和 MNTD 标准化效应量的相反数：`NRI = (零模型 MPD 均值 - 观测 MPD) / 零模型 MPD 标准差`，NTI 对应使用 MNTD。正值表示相对于所选零模型的系统发育聚集，负值表示离散。

对于单个观测群落及其本地物种池，使用 `nri_nti`：

```julia
pool = ["A", "B", "C", "D"]
observed = ["A", "B", "D"]

result = nri_nti(tree, pool, observed; nreps = 999, seed = 20260713)
result.nri
result.nti
result.status   # :ok 或 :null_sd_zero
```

零模型从物种池中无放回抽取与观测群落物种数相同的物种。支持 tip 标签和 tip 位置（`1:tree.ntips`）。`nri_nti` 使用稳定随机数生成器和局部虚拟树计算 MPD/MNTD，不假设全球物种池，也不在内存中保存全部零模型重复结果。兼容属性还包括 `NRI`、`NTI`、`MPD_obs` 和 `MNTD_obs`。

批量分析空间网格时，可传入 `grid_id`，使每个网格的随机种子由全局种子稳定派生，不受处理顺序影响：

```julia
result = nri_nti(tree, pool, observed;
    nreps = 999,
    seed = 20260713,
    grid_id = "1745",
)
```

派生后的实际种子保存在 `result.seed`。不要同时传入 `grid_id` 和显式的 `rng`，因为传入的 RNG 已定义随机数流。

对于大量网格，`nri_nti_table` 接受两张标准分隔长表：一张为本地物种池，一张为观测群落；两表均使用 `grid_id`、`species` 两列。表格可以是 CSV、TSV、TXT 或 CSV.jl 支持的其他分隔格式；需要显式指定制表符时使用 `delim = '\t'`。每个文件中同一网格的记录必须连续。Julia 启动了多个线程时可用 `threaded=true` 并行处理网格。下例临时创建一个网格的输入文件：

```julia
batch = mktempdir() do dir
    pool_path = joinpath(dir, "local_pool.tsv")
    observed_path = joinpath(dir, "observed_community.tsv")
    write(pool_path, "grid_id\tspecies\ngrid_1\tA\ngrid_1\tB\ngrid_1\tC\ngrid_1\tD\n")
    write(observed_path, "grid_id\tspecies\ngrid_1\tA\ngrid_1\tB\ngrid_1\tD\n")
    nri_nti_table(tree, pool_path, observed_path;
        nreps = 99, seed = 20260713, threaded = false, delim = '\t')
end
batch.grid_ids
batch.results
```

## 系统发育 beta SES

`phylobeta_ses` 复用现有 SES 的稀疏随机化和在线统计，对一个系统发育 beta 成分计算标准化效应量。每次选择 `:sorensen` 或 `:jaccard`，以及 `:total`、`:turnover` 或 `:nestedness` 之一。

```julia
using Random

effect = phylobeta_ses(data, tree; index_family = :sorensen,
    component = :turnover, null_model = :taxa_labels, runs = 999,
    pairs = [(1, 2)], rng = MersenneTwister(42))
effect.observed
effect.null_mean
effect.null_sd
effect.z_score
effect.p_value
effect.valid_replicates
```

支持稠密/CSC 矩阵和 `CommunityData`，所有正丰度视为存在。完整树物种池始终参与抽样，包括未观测物种；`CommunityData` 按物种名对齐并补零，直接矩阵输入须保留每个树 tip 的列。

随机化遵循[共用零模型表](#选择零模型)。

`PhylobetaSESResult` 包含上述六组数组、群落标签和丰富度、指定对及零模型参数。默认返回完整对称矩阵；`pairs` 返回按请求顺序排列的向量，保留重复、反向及自配对。每次仍随机化完整输入，仅构建指定端点的分支集合；在线统计不保留各次随机距离矩阵，指定对输出也不分配完整群落间矩阵。

均值和样本 SD 只使用有限随机值，`valid_replicates` 逐对记录次数。没有有效值时均值/SD/z/p 为 `NaN`；只有一个时 SD/z 为 `NaN`；常量零模型保留 SD=0、z=`NaN`。观测 beta 非有限时，即使有有限的零模型统计，z/p 仍为 `NaN`。空群落和零 PD 群落沿用 `phylobeta` 规则，包括未定义的对角线。

下尾中秩分数为 `p = (1 + 小于观测值的次数 + 0.5 × 相等次数) / (有效次数 + 1)`，相等按计算值精确判定。它不是双尾 p 值；分母采用有效次数，群落内 SES 使用 `runs + 1`。R 的 `rowwise`、`colwise` 标签重排会被对齐抵消，因此拒绝这两个模型；也不沿用 R 将零 SD 替换为一的做法。SES z 值可能为负或未定义，不宜作为 `phyloregion` 的距离输入。
