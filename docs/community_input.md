# Community inputs and alignment

[Documentation](README.md) · [Home](../README.md) · [简体中文](community_input_zh.md)

Keep community values and labels together with `CommunityData`. Rows are sites,
columns are species, and values must be finite nonnegative real numbers.
Positive values denote presence; abundance-weighted functions use their magnitudes.

## Matrices, long tables and triplets

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

| Input | Storage and labels |
| --- | --- |
| Matrix | Copies values and retains dense/CSC storage and numeric type; `storage=:sparse` requests CSC. |
| Long `DataFrame` | Builds CSC directly; labels default to first occurrence. `value=nothing` counts records. |
| Indexed triplets | Builds CSC directly; provide complete `sample_labels` and `species_labels`. |

Explicit label lists preserve empty sites, all-zero species columns and the
desired order. Labels must be unique nonempty strings or integer identifiers;
integers become strings. Matrix defaults are `sample_1`, … and `species_1`, …;
provide actual species names for tree-based work. Missing or nonfinite values fail.

Long-table and triplet duplicates sum by default. `duplicates=:error` rejects
them, including zero-valued duplicates; overflowing sums fail. Integer abundances
remain counts, fractional weights remain fractional, and Bool records become Int.

## Export a long table and reconstruct the community

`community_longtable` exports one row per positive cell, in species-then-site
order. Keep its two complete label lists: positive records alone lose empty
sites and all-zero columns. Values are not expanded into repeated occurrences.

```julia
exported = community_longtable(data)
exported.table                     # sample, species, value columns
restored = CommunityData(exported.table;
    sample_labels = exported.sample_labels, species_labels = exported.species_labels)
restored.matrix == data.matrix      # true
restored.sample_labels == data.sample_labels   # true
restored.species_labels == data.species_labels # true
```

This also preserves zero-sized domains and entirely empty tables. Exported
arrays are independent; reconstruction builds CSC and follows the same numeric
rules as other long tables. See [export costs](performance.md#community-long-table-export).

## Tree alignment and metric calls

```julia
tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
aligned = align_community(data, tree)
aligned.species_labels             # ["A", "B", "C", "D"]
pd(data, tree).pd                   # label alignment is automatic
phylo_endemism(data, tree).pe
```

`CommunityData` fills absent tree species with zero columns and preserves the
tree and root. Unknown input species always fail, even in an all-zero column.
Set `missing_tips=:error` in `align_community` to require every tree species.
Alignment preserves dense/CSC storage. Set labels at construction; do not also
pass `tip_labels`/`sample_labels` to metric calls on `CommunityData`.

Direct matrix calls require one column per tree tip, including unobserved tips.
Without `tip_labels`, columns use `tree.tip_labels`; otherwise the keyword
describes the supplied column order. Results retain input site order.

## Repeated analyses

| Repeated work | Prepare once | What is reused |
| --- | --- | --- |
| Species matching only | `align_community(data, tree)` | Aligned values and labels |
| PD, WE/PE, RPD, UniFrac/PhyloSor, beta | `prepare_community(data, tree)` | Branch sets, occupancy and per-site summaries |
| `comdist` / `comdistnt`, either weighting mode | `prepare_community_distances(data, tree)` | Abundance/presence records, working topology and distance workspaces |
| Cuts of the same clustering hierarchy | `prepare_phyloregion(beta)` | Distance snapshot and hierarchy |

```julia
ctx = prepare_community(data, tree)
pd(ctx).pd
phylo_endemism(ctx).pe
phylobeta(ctx; pairs = [(1, 2)])
```

Constructors, alignment and contexts own their arrays. Treat data/context fields
as read-only; changing original inputs does not update a snapshot. Prepare again
for changed inputs. Metric results can be modified independently. Contexts retain
different summaries, so use the context supported by the function being called;
they do not cache every result or every full pairwise matrix.

Continue with [diversity and beta](community.md), [community distances](community_distances.md)
or [null models](null_models.md). [API conventions](api.md) summarize result shapes.
