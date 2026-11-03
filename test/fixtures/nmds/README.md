# NMDS reference configurations

`r_reference.tsv` is one compact table generated offline on 2026-11-03.
Four rows cover a tied taxonomic beta input and duplicated sites, each with
weak and strong ties. Every row stores the complete distance matrix, starting
coordinates, fitted Euclidean distances, monotone disparities and stress-1.
Pair vectors use upper-triangle column order: (1,2), (1,3), (2,3), ... .

The unmodified local vegan 2.8-0 `monoMDS.R` wrapper was evaluated with
`maxit=0`, scaling enabled and PC rotation disabled, using the installed vegan
2.7-5 Fortran engine in R 4.6.2. Results were identical to the installed wrapper
for these calls. These fixtures compare the objective at the same configuration;
SMACOF and vegan's steepest-descent optimizer need not follow the same path.
The taxonomic beta validation input is explicitly rounded to 12 decimal places
to avoid vegan's absolute near-tie tolerance differing from Julia's exact ties.
Production `nmds` never rounds inputs.

Ordinary Julia tests only read this file and require neither R nor internet.
The larger eight-case comparison and logs remain under ignored
`validation/community/`; reproduction uses `nmds_export.jl`,
`nmds_reference.R`, `nmds_crosscheck.jl` and `nmds_check.jl` there.
