# Phylogenetic beta SES reference data

`r_reference.tsv` is one small fixed table generated offline on 2026-11-03
with local R 4.6.2, phyloregion 1.0.9 source, Matrix, ape, phangorn and betapart.
Ordinary Julia tests read this file and never start R or install R packages.

Four examples retain the tree, complete community, null parameters and
StableRNG seed with R reference summaries. They cover taxa-label shuffling,
frequency randomization with empty sites, a zero-length tree and a fixed-margin
swap null. Expected output includes every lower-triangle entry and its finite
replicate count; undefined-site diagonals retain NaN.

The R generator evaluates independently drawn dense null communities using the
real local phyloregion `phylobeta` / betapart, then computes means and sample
SDs with base R. It does not translate Julia RNG seeds to R RNG seeds: both
consume the same exported reference communities. Exact last-bit comparisons
can classify rank ties differently across languages; `p_min`/`p_max` bound this
using 32 machine eps, while the Julia dense-draw oracle checks exact ranks.

The full local comparison covers six data sets, seven null models, two beta
families and three components (252 combinations), plus three calls of the
unmodified original `phylobeta_ses` wrapper with aligned tip permutations.
Its rowwise and colwise degeneracy and zero-SD replacement are separately
validated. Scripts live in the ignored `validation/community/` directory:
`beta_ses_export.jl`, `beta_ses_reference.R`, `beta_ses_crosscheck.jl`.
