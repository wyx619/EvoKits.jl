# Conservation metrics

[Home](../README.md) · [简体中文](conservation_zh.md)

EDGE2 scores and species extinction probability tools.

The example below uses this four-tip tree:

```julia
using EvoKits

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
```

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
