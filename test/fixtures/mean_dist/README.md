# Mean distance reference data

`r_reference.tsv` contains four small fixed ensembles generated offline on
2026-11-03 with R 4.6.2, Matrix 1.7-6, ape, phangorn, betapart 1.6.1 and local
phyloregion 1.0.9. Ordinary Julia tests read this one file without starting R.

Cases cover finite matrices with reordered labels, subsets, NaN observations,
and one posterior-tree phylogenetic beta ensemble. Each lower-triangle input
row includes its local row indices/labels and the output label order. Expected
means/counts are repeated so inputs and reference results stay in one table.

`mean_skip`, `mean_propagate` and `count` are independently computed in R from
the observed values for each labeled pair. `native_mean` comes from the
unmodified local `mean_matrix.R` function on R `dist` inputs: it propagates
NaN, sorts the first input's labels, and materializes diagonals as zero. Its
off-diagonal values are compared under matching propagation rules. Julia
preserves undefined diagonals and reports finite counts explicitly.

The larger local validation covers 42 ensembles, including 18 real posterior
beta ensembles (three data sets, two index families, three components). Its
generator and checker are ignored local files:
`validation/community/mean_dist_reference.R` and `mean_dist_crosscheck.jl`.
Only these four examples are kept in the ordinary test fixture.
