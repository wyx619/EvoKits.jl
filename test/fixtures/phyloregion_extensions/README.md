# Small phyloregion reference

`reference.tsv` contains 24 fixed cases generated on 2026-11-03 with R 4.6.2
from the local phyloregion 1.0.9 sources: `functional_beta.R` (`lbow`),
`Choropleth.R` (`JENKS`) and type-7 upper quantiles used by `coldspots.R`.
It stores four elbow curves, eight Jenks partitions and twelve quantile cases.
All inputs are embedded in the table. Normal Julia tests use no R process,
network connection or optional R packages.

Jenks partitions can tie: tests compare the minimum within-class SSE rather
than requiring identical breakpoints. Julia keeps duplicate values together,
retains zeros and fractions, and uses scaled Welford interval variances.
Hotspot threshold tests compare quantiles directly. The R coldspot function's
zero-threshold mutation error is not adopted.

The ignored local reproduction scripts are
`validation/community/extensions_reference.R` and `extensions_crosscheck.jl`.
Their larger comparison covers 437 outputs, including KL rankings, mixed-trait
group-count selection, EDGE and tree placement distances. Only this small
table is needed by `test/phyloregion_extensions.jl`; independent formulas and
boundary cases cover the other features.
