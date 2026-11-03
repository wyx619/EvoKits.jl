# Infomap reference

One TSV contains 18 site rows from three small, stable binary networks: two
disconnected blocks, three disconnected blocks, and a complete bipartite graph.
Each row includes its input values and the reference site assignment.

Generated from the local phyloregion 1.0.9 `infomap` wrapper with R 4.6.2,
R igraph 2.3.3 / C igraph 0.10.17, seed 42 and 10 trials. Partitions are
renumbered by the first input site. Ordinary tests read this table offline;
R and the ignored validation/reference directories are unnecessary.

The R core has no regularization option and records teleportation in its flow
model even for undirected graphs. Tests compare these stable site partitions
with `regularized=false`. The R code lengths are provenance, not expected
scores for the newer engine or evidence that default regularization matches R.
