# Multiple-site taxonomic beta reference

`reference.tsv` contains eight outputs: four small presence/absence collections
under the Sørensen and Jaccard families. Matrices are packed in column-major
order; `NaN` records undefined components.

Generated locally with R 4.6.2, betapart 1.6.1 and Matrix 1.7-6. Each output
agrees between `betapart::beta.multi` applied to the unmodified phyloregion
1.0.9 `beta_core` and `betapart::beta.multi` applied directly to the matrix.
Cases cover mixed turnover/nestedness, strict nesting, an empty/nonempty pair,
and all-empty sites. Julia separately rejects fewer than two selected sites.

The generator, larger comparisons and R session information stay under ignored
`validation/`; ordinary tests only read this single table and do not start R.
