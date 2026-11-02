# Functional grouping reference data

These fixed test inputs and reference outputs were generated offline on
2026-11-03 with R 4.6.2, clustMixType 0.5.2, Matrix 1.7-6 and betapart 1.6.1,
using local phyloregion 1.0.9 source. Ordinary Julia tests read the files and
do not start R or require any R packages.

Only six representative cases from the 40-case comparison are retained here:
one hand calculation, a single group, all three numeric scales, unequal trait
weights, explicit/automatic gamma, and a constant numeric column.

- `kproto_r.tsv` contains input traits, settings and R membership, objective,
  numeric prototypes and categorical modes. R uses exactly the initial species
  selected by the saved Julia StableRNG seed, explicit weights, and independently
  transformed numeric traits. `gamma_mode=auto` is independently resolved in R
  from all species pairs; the R package's default lambda estimate differs.
- `community_r.tsv` contains integer input abundances and R sparse aggregation
  values. Outputs repeat per input species so the complete input and mapping
  remain in one table. Presence, richness, abundance and abundance proportions
  include an empty site; its proportions are explicitly set to Julia's zero
  policy after R's otherwise undefined division.
- `beta_r.tsv` contains lower-triangle Sørensen/Jaccard components from local
  phyloregion `beta_diss` and betapart, including empty-site comparisons. Empty
  diagonal entries are not compared because R `dist` materializes them as zero
  while Julia reports NaN.

The ignored scripts `validation/community/functional_export.jl`,
`functional_reference.R` and `functional_crosscheck.jl` reproduce all 40 cases.
The larger local comparison also executes the unmodified R `functional_beta`
body with injected calls controlling only k-prototypes starts/weights and
capturing its aggregation, plus a complete call with no injected dependencies.
Random-number sequences, tie assignment, duplicate/empty-group handling and
automatic gamma defaults are not claimed to match between languages.
