# Conservation metrics

[Documentation](README.md) · [Home](../README.md) · [简体中文](conservation_zh.md)

On this page: [EDGE](#edge) · [EDGE2](#edge2) · [Hotspots/coldspots](#hotspots-and-coldspots) · [Value classes](#metric-value-classes)

EDGE/EDGE2 scores, extinction probabilities, metric tail selection and value classification.

The example below uses this four-tip tree:

```julia
using EvoKits

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
```

## EDGE

`edge` combines full-tree fair-proportion ED and a category-derived GE score:
`log1p(ED) + GE*log(2)`. Its default mapping follows phyloregion's actual
probability-based implementation:

| Category | LC | NT | VU | EN | CR | EW/EX | DD |
|---|---|---|---|---|---|---|---|
| GE | 0.001 | 0.01 | 0.1 | 0.67 | 0.999 | 1 | missing |

```julia
edge_result = edge(tree, ["LC", "NT", "VU", "CR"])
edge_result.species_labels
edge_result.ed
edge_result.edge
edge_result.status

# Choose integer GE category scores explicitly when that convention is needed.
integer_edge = edge(tree, ["LC", "NT", "VU", "CR"];
    ge_mapping = Dict("LC" => 0, "NT" => 1, "VU" => 2, "EN" => 3, "CR" => 4))
```

The probability mapping differs from the common integer category coding;
compare results only under the same mapping. A custom dictionary replaces the
default and must cover every supplied category. `EDGEResult` returns full-tree
tip order, ED, GE, EDGE, status and the mapping used. A vector uses tree order
unless `tip_labels` supplies exact species labels; a `DataFrame` uses
`species`/`category` columns. Subsets are accepted. Absent records, missing
categories and DD return missing scores by default, with `:absent_record`,
`:missing_category` or `:unassessed`; `missing_category=:error` rejects them.
Unknown species/categories and duplicate labels always fail. Branch lengths
must be finite and nonnegative. EDGE2 below has its own model and result.

## EDGE2

EDGE2 combines evolutionary distinctiveness with species extinction
probabilities. It uses the same `CompactTree` representation as the other
EvoKits tools, so no separate tree parser or `Phylo.jl` dependency is needed.

```julia
pext = [0.10, 0.20, 0.30, 0.40]  # tree.tip_labels order
edge2 = edge2_mod(tree, pext)

edge2.species
edge2.edge       # EDGE2 scores
edge2.ed         # evolutionary distinctiveness
edge2.pd         # total phylogenetic diversity
edge2.epd_loss   # expected PD loss
```

The reference field spellings `EDGE`, `ED`, `TBL`, `PD`, and `ePDloss` remain
available as compatibility properties. Extinction probabilities must be
finite values in `(0, 1]` and must follow the tree tip order.

For category-based workflows, use `PextDistribution`,
`calc_pext_distribution`, `load_pext_distribution_csv`, and `sample_pext`.
`reorder_pext_to_tree` validates species coverage before returning values in
tree order.

## Hotspots and coldspots

`hotspots` selects values at or above an upper quantile; `coldspots` selects
values at or below a lower quantile. `prob` is a tail percentage, default 2.5.
Both use R type-7 interpolation, include all threshold ties, and retain input
order and labels. The selected fraction can exceed the requested percentage.

```julia
metric = [0.0, 1.0, 2.0, 10.0, 12.0]
high = hotspots(metric; prob = 20, sample_labels = ["a", "b", "c", "d", "e"])
low = coldspots(metric; prob = 20)
high.threshold
high.selected

# A common reference keeps thresholds comparable between datasets.
shared_threshold = hotspots([1.0, 8.0, NaN]; reference = metric,
    prob = 20, invalid = :skip)
shared_threshold.valid_count
shared_threshold.reference_count
```

`hotspots(pd_result)`, `hotspots(pe_result)` and `hotspots(we_result)` also
accept the corresponding metric results and preserve site labels; so do
`coldspots`. `MetricSelectionResult` includes values, flags, threshold, tail,
percentage and valid counts. `invalid=:skip` preserves missing/NaN positions
with missing flags; Inf always fails. The threshold reference needs a valid
value. Zero and negative metric values participate normally. Selection applies
the quantile inequality directly, correcting R `coldspots`' zero-threshold
mutation error.

## Metric value classes

`classify_values` provides exact Jenks natural breaks (`:jenks`, default),
equal widths (`:equal`) and type-7 quantiles (`:quantile`).

```julia
classes = classify_values(metric; k = 3, method = :jenks)
classes.classes       # original input order; classes start at 1
classes.upper_bounds
classes.breaks
classes.objective     # Jenks within-class SSE
classes.actual_k
```

`ValueClassificationResult` also retains labels, values, requested count,
method and reason. Classes include their upper boundary; repeated values
always share a class. `breaks` contains the minimum followed by all upper
bounds, so a singleton minimum class gives two equal initial entries. Jenks
uses at most the distinct-value count; equal/quantile methods cap bins at the
valid observation count and deduplicate upper bounds. `actual_k` counts those
bounds; equal/quantile intervals can be empty. A count below k gives
`:reduced_classes`; constants give one class. Empty/all-invalid inputs fail.
`invalid=:skip` leaves missing classes at missing/NaN positions; Inf fails.

Jenks minimizes total within-class SSE, accounting for repeated-value counts.
`max_work=50_000_000` bounds the exact workload; excess work fails explicitly.
Original-unit SSE may overflow/underflow; `normalized_objective`
and `objective_scale` retain the scaled result. Equal-cost splits prefer the
latest final split recursively, which may differ from R while giving the same
minimum SSE. Values retain zeros, negatives and fractions. These are numeric
classes; map geometry, colors and plotting are separate operations. See
[exact-classification costs](performance.md#metric-value-classification).
