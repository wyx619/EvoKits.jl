# Phylogenetic multiple-site beta reference

`reference.tsv` contains eight results for four small collections and both
Sørensen/Jaccard families. It was generated offline with R 4.6.2, betapart 1.6.1
and the local phyloregion `R/phylo_beta.R` on 2026-11-04.

The phyloregion core's `St` was replaced by the PD of branches occupied by the
selected sites. Each result was also checked against betapart's unmodified
`phylo.betapart.core` + `phylo.beta.multi` tree calculation. The unobserved case
retains the original root path while excluding unoccupied branches; empty
cases record undefined components as `NaN`.

Matrices are column-major comma-separated presence values; `columns` gives
species order and `selected` the original one-based row indices. `newick`
contains the tree. The three components are scalars, and `pooled_pd` is the
selected branch-union length. No R installation is needed for normal tests.

The full generator, original-core differences, randomized comparisons and R
session details remain in the ignored local `validation/community/` directory.
