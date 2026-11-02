# Tree lineage reference data

These are fixed regression inputs and outputs, not analysis results written by
normal package use. Ordinary Julia tests read the TSV files without starting R.

- `trees.tsv`: 15 Newick trees, including integer-age, multifurcating,
  zero-terminal-edge, and random coalescent examples.
- `partitions_r.tsv`: full tip membership from local phyloregion `get_clades.R`,
  for interior cuts, selected exact boundaries and all requested k values.
- `communities_r.tsv`: input abundances and expected lineage presence from
  local `collapse_range.R`; richness/abundance from the sparse multiplication
  core of `counts.R`. Output values are repeated per input tip to keep the
  original species columns explicit.

Generated offline on 2026-11-02 with R 4.6.2, ape, phangorn and Matrix, using
seed 20261106. No installation or network calls. Generator and session details
are kept locally in ignored `validation/community/lineage_reference.R` and
`results/lineage_r_session.txt`. Reference source is local phyloregion 1.0.9.

R comparisons cover shared semantics. Julia separately tests the root and
all-zero boundaries, explicit younger-side cuts, full-tree alignment and
preservation of empty sites/groups. R's pruning, automatic name changes,
fixed k epsilon and long-table conversion are not adopted.
