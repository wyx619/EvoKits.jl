# Community phylogenetic distances

[Documentation](README.md) · [Home](../README.md) · [简体中文](community_distances_zh.md)

Compare communities by average or nearest-taxon distance, or use branch overlap.
Construct and align inputs as described in [community inputs](community_input.md).

```julia
using EvoKits, SparseArrays

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end
community = [1 0 1 0; 0 1 1 0]
data = CommunityData(sparse(community);
    sample_labels = ["west", "east"], species_labels = tree.tip_labels)
```

## Distances between communities

`comdist` averages all cross-community tip distances. Its diagonal includes
zero self-distances as well as distances between different tips. Positive
values supply equal species weights by default; `abundance_weighted=true`
uses their abundances. Empty sites give `NaN` rows and columns.

`comdistnt` averages the nearest-target distance in each direction. Directions
receive equal weight by default; with `abundance_weighted=true`, source
abundances weight each directional mean and original site totals weight the
two means. A source tip without an eligible target is omitted and remaining
weights are renormalized. Use the valid direction if only one exists, or
`NaN` if neither does. Nonempty diagonals are zero by default.

```julia
distance_data = CommunityData(sparse([2.0 1 0; 0 1 3; 0 0 0]);
    sample_labels=["west", "east", "empty"], species_labels=["A", "B", "C"])
cross_mean = comdist(distance_data, tree; abundance_weighted=true)
nearest = comdistnt(distance_data, tree; exclude_conspecifics=true)
cross_mean[3, :] # [NaN, NaN, NaN]
nearest[1, 1]   # 2.0: nearest different species within the first site
requested_pairs = [(1, 2), (2, 1), (1, 1), (3, 2)]
selected_mean = comdist(distance_data, tree;
    abundance_weighted=true, pairs=requested_pairs)
# [10/3, 10/3, 8/9, NaN]
selected_nearest = comdistnt(distance_data, tree;
    abundance_weighted=true, exclude_conspecifics=true, pairs=requested_pairs)
# [22/7, 22/7, 2, NaN]
```

Both `comdist` and `comdistnt` accept `pairs` and return a vector in request order, including
repeated, reversed and self pairs. Pair indices refer to the complete input
site order; an empty request returns an empty vector. Values, labels and tree
lengths are validated even for sites absent from the requests. Only requested
endpoints are retained for computation, and repeated/reversed requests reuse
one evaluation. `comdistnt` computes both directional means for each unique
pair, reusing each needed target pass. To average off-diagonal pairs with `mean_dist`, supply labeled
tuples containing `distances`, `sample_labels` and `pairs`, and explicit
output `pairs`. `mean_dist` requires self-distances to be zero or NaN, while
`comdist` self-pairs include within-site diversity. `comdistnt` self-pairs are
zero by default and may be nonzero with `exclude_conspecifics=true`.

`exclude_conspecifics=true` excludes species identity. Different species can
still have distance zero on zero-length branches; singleton diagonals then
give `NaN`. These are the existing Julia conventions. Picante combines
directional means differently when richness or site totals differ and excludes
all zero distances; it also returns a `dist` object with zero diagonals.

Both functions accept dense/CSC matrices or `CommunityData`, preserve site
order, and require finite non-negative branch lengths with a finite sum.
They leave inputs unchanged and support column alignment through `tip_labels`
or `CommunityData` labels. Their CSC paths avoid a dense community and a full
tip-distance matrix. Default output is a full site-pair matrix; selected
distance output contains only the requested values. See
[algorithm and storage costs](performance.md#community-phylogenetic-distances).

For repeated distance calculations, prepare a shared snapshot:

```julia
distance_context = prepare_community_distances(distance_data, tree)
reused_mean = comdist(distance_context; abundance_weighted=true,
    pairs=requested_pairs)
reused_nearest = comdistnt(distance_context; abundance_weighted=true,
    exclude_conspecifics=true, pairs=requested_pairs)

fixed_context = prepare_community_distances(distance_data, tree;
    pairs=requested_pairs)
fixed_nearest = comdistnt(fixed_context; abundance_weighted=true,
    exclude_conspecifics=true) # same request order as selected_nearest
```

`CommunityDistanceContext` stores both abundance and presence weights, one
owned tree topology and the workspaces for both distances. With no preparation
`pairs`, later calls can request full matrices or change the pairs; these calls
reuse the topology for all prepared sites. Fixing `pairs` during preparation
retains only their nonempty endpoints and caches request schedules. Such calls
default to the fixed vector request; rebuild to change its scope or order.

All input values, species labels, branch lengths and site labels are checked
at preparation; site labels must be unique and nonempty. The original matrix
and tree are not retained, and later edits do not affect the snapshot. Rebuild
for changed data. Treat context fields as internal and read-only. Results own
their arrays; calls sharing a context are serialized to protect its reusable
workspaces. No distance matrix is cached. The presence-based `CommunityContext`
for PD/PE/beta retains different summaries and remains a separate entry point.


## UniFrac and PhyloSor

`unifrac(data, tree)` returns unweighted branch-length dissimilarity;
`phylosor(data, tree)` returns similarity `2*shared_PD/(PD1+PD2)`.
Both use presence, preserve row order and return complete matrices.

UniFrac returns zero when the combined branch length is zero; a zero-PD site
versus a positive-PD site gives one. PhyloSor gives one on a positive-PD site's
diagonal and zero whenever its denominator is zero, including empty and
zero-PD sites. See [matrix semantics](api.md#matrix-semantics)
for zero-denominator behavior. Use `phylobeta` when turnover/nestedness
decomposition or selected branch-overlap pairs are needed.

```julia
unifrac(data, tree)
phylosor(data, tree)
phylobeta(data, tree; pairs = [(1, 2)])
```

Clustering requires a complete distance matrix with finite values and zero
diagonal. `comdist` includes within-site diversity on the diagonal and PhyloSor
is a similarity; interpret and transform these explicitly before using them
as distances. [Regionalization](regionalization.md) documents accepted inputs.
