# NMDS postprocessing reference

`reference.tsv` holds four small configurations evaluated locally with R 4.6.2:
two labeled Procrustes transforms (including a reflection and uniform scaling)
and two coordinate-to-color mappings with different aspect ratios.
Matrices are packed in column-major order; coordinates use 17 significant digits.

Alignment uses the unmodified `procrustes.R` from the local vegan 2.8-0 source,
with vegan 2.7-5 providing score extraction. Colors use the unmodified local
phyloregion 1.0.9 `hexcols_rescale.R` and colorspace 2.1-3, with numeric R row
names; arbitrary sample labels are preserved separately by EvoKits.

Generators, additional fits, all 40,000 palette cells and session information
remain under ignored `validation/`. Ordinary Julia tests only read this table.
