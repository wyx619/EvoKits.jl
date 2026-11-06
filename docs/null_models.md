# Null models and community structure

[Documentation](README.md) · [Home](../README.md) · [简体中文](null_models_zh.md)

On this page: [Site SES](#site-ses) · [Null model choice](#choose-a-null-model) · [NRI/NTI](#nri-and-nti) · [Beta SES](#phylogenetic-beta-ses)

SES and NRI/NTI compare community structure with a chosen null model. They
describe species assemblages; trait-based phylogenetic signal is covered
in [phylogenetic signals](signals.md).

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

## Site SES

```julia
using Random

pd_effect = ses_pd(data, tree; null_model = :richness,
    runs = 99, rng = MersenneTwister(42))
mpd_effect = ses_mpd(data, tree; null_model = :richness,
    runs = 99, rng = MersenneTwister(42))
mntd_effect = ses_mntd(data, tree; null_model = :richness,
    runs = 99, rng = MersenneTwister(42))
pd_effect.pd_obs_z
mpd_effect.mpd_nri             # negative of mpd_obs_z
mntd_effect.mntd_nti           # negative of mntd_obs_z
```

Each result retains input site order, richness (`ntaxa`), observed values,
null means and sample SDs, ranks, z scores and p-values, plus run settings.
Field prefixes are `pd_`, `mpd_` and `mntd_`. Use `ses_pd_dataframe`,
`ses_mpd_dataframe` or `ses_mntd_dataframe` for tabular output.

`z = (observed - null_mean)/null_sd`. Empty MPD/MNTD observations are NaN;
singletons are NaN except weighted MPD with positive mass, which is zero.
PD follows the [root/singleton rules](api.md#boundary-policies).
Abundance-weighted MPD/MNTD randomize presence while keeping original abundances
at their tree-tip positions; zero-weight null tips remain MNTD candidates.

Means/SDs use finite null draws. With none, mean/SD/z are NaN; with one, SD/z
are NaN; zero SD also gives NaN z. The lower-tail mid-rank score is
`(1 + count(null < observed) + 0.5*count(null == observed))/(runs + 1)`.
Nonfinite draws add no rank count. It is not a two-sided p-value.

## Choose a null model

| Null model | Randomization and constraints |
| --- | --- |
| `:taxa_labels` (`:tipshuffle`) | One shared permutation of tree-tip labels; preserves taxonomic overlaps. |
| `:richness`, `:sample_pool` | Independently samples each site's richness from the complete tree-tip pool. |
| `:phylogeny_pool` | Richness sampling followed by a shared label permutation. |
| `:frequency` | Preserves each species' site occupancy; site richness can change. |
| `:trialswap`, `:independentswap` | Checkerboard trials preserve both margins; each replicate restarts from the observed community. |

The complete tree-tip pool includes unobserved species. The `:sample_pool`
alias uses the same pool as `:richness`. Swaps use `iterations` and
`max_attempts_per_swap`; finite swap effort does not guarantee uniform draws.
Use fresh RNGs with the same seed/settings to compare statistics on matching
draws. Null draws are summarized online; replicate communities are not retained.

## NRI and NTI

NRI and NTI describe within-community phylogenetic structure relative to a
null model. They are the negative standardized effects of MPD and MNTD:
`NRI = (null_mpd_mean - observed_mpd) / null_mpd_sd`, with the corresponding
MNTD formula for NTI. Positive values indicate clustering and negative values
indicate overdispersion relative to the chosen null model.

For a single observed community and a local species pool, use `nri_nti`:

```julia
pool = ["A", "B", "C", "D"]
observed = ["A", "B", "D"]

result = nri_nti(tree, pool, observed; nreps = 999, seed = 20260713)
result.nri
result.nti
result.status   # :ok or :null_sd_zero
```

The pool is sampled without replacement at the observed richness. Both tip
labels and tip positions (`1:tree.ntips`) are accepted. `nri_nti` uses a stable
random generator and a local virtual-tree MPD/MNTD calculation; it does not
assume a global species pool or retain all null replicates in memory. The
reference field names `NRI`, `NTI`, `MPD_obs`, and `MNTD_obs` are also available
as compatibility properties.

For a batch of spatial grids, pass `grid_id` so the global seed is stably
derived per grid, independent of processing order:

```julia
result = nri_nti(tree, pool, observed;
    nreps = 999,
    seed = 20260713,
    grid_id = "1745",
)
```

The effective derived seed is stored in `result.seed`. Do not combine `grid_id`
with an explicit `rng`; an injected RNG already defines the random stream.

For many grids, `nri_nti_table` accepts two standard delimited long tables with
the columns `grid_id,species`: one local-pool table and one
observed-community table. Tables may be CSV, TSV, TXT, or another delimiter
recognized by CSV.jl; use `delim = '\t'` when an explicit TSV delimiter is
needed. Each grid must be a contiguous block in each file. Grids are processed
in parallel with `threaded=true` when Julia starts with multiple threads.
This example creates one grid's input files temporarily:

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

## Phylogenetic beta SES

`phylobeta_ses` standardizes one phylogenetic beta component against a null
model. It reuses the sparse randomization and online statistics of the site
SES methods. Choose `:sorensen` or `:jaccard` and one of `:total`, `:turnover`,
or `:nestedness` per call.

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

Dense/CSC matrices and `CommunityData` are accepted; positive abundances
denote presence. The complete tree-tip pool remains eligible, including
unobserved species. `CommunityData` aligns by species name and fills missing
tree species with zeros; matrix inputs require a column for every tree tip.

Randomization follows the [shared null-model table](#choose-a-null-model).

`PhylobetaSESResult` contains the six arrays shown above, site labels and
richness, requested pairs, and null-model settings. Arrays are symmetric full
matrices by default. `pairs` returns vectors in request order, retaining
repeated, reversed, and self pairs. Every draw randomizes the complete input;
only branch sets for requested endpoints are built. Replicate matrices are
never retained, and selected-pair output avoids full site-by-site matrices.

Means and sample SDs use only finite draws, counted per pair. No valid draw
gives `NaN` mean/SD/z/p; one valid draw gives `NaN` SD/z. Constant nulls retain
SD zero and yield `NaN` z. Nonfinite observed beta gives `NaN` z/p even when
finite null summaries exist. Empty and zero-PD sites follow `phylobeta`,
including undefined diagonals.

The lower-tail mid-rank score is
`p = (1 + count(null < observed) + 0.5*count(null == observed)) / (valid + 1)`.
Ties use exact computed values. This is not a two-sided p-value; its
denominator uses valid draws, while site SES uses `runs + 1`.
R's `rowwise` and `colwise` label reorderings are undone by alignment, so
these models are rejected. R's replacement of zero SD by one is not used.
SES z scores can be negative or undefined and are unsuitable as distances
for `phyloregion`.
