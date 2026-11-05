# Trees

[Home](../README.md) · [简体中文](trees_zh.md)

Read, query, prune, complete and map trees, or compare branch history at different time depths.

## Tree I/O and Topology

```julia
using EvoKits, SparseArrays

tree = mktempdir() do dir
    path = joinpath(dir, "example.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    read_tree(path)
end

ultrametric_tree = make_ultrametric(tree)
binary_tree = multi2di(ultrametric_tree; random = false)
write_tree("binary.nwk", binary_tree)
```

`CompactTree` is the common input type for EvoKits tools. It caches topology,
edge lengths, traversal order, and tip-label lookup. `read_tree` reads a Newick
file directly into this representation.

Tree queries use compact node IDs, tip positions, or tip labels explicitly:

```julia
mrca(tree, "A", "B")
distance(tree, 1, 2)
pairwise_distance(tree)
descendant_tips(tree, tree.root; labels = true)
tree_imbalance(tree; metric = :sackin)
```

Weighted distances use branch lengths; pass `weighted = false` for edge
counts. `distance`, `pairwise_distance` and `distances_to_tip` preserve short
paths below long shared stems by summing original branch segments when
root-depth subtraction would lose precision. Tree code is grouped into
queries, operations, mapping and metrics;
see the [source layout](development.md) for file locations and dependency order.

## Posterior tree streams and samples

`read_trees` reads Newick or NEXUS files one tree at a time, detecting each
file's format from its contents. Newick records end
at semicolons and may span lines or share a line; bracket comments, rooting
annotations and single-quoted labels are supported. Supply one path or a
vector/tuple of paths in the desired order. `burnin` discards initial records
and `thin` retains every nth remaining tree, independently within each file.

```julia
using Random

posterior_example = mktempdir() do dir
    path = joinpath(dir, "posterior.nwk")
    write_tree(path, "((A:1,B:1):1,(C:1,D:1):1);")
    write_tree(path, "((C:1,D:1):2,(B:1,A:1):2);"; append = true)
    write_tree(path, "((A:1,C:1):1,(B:1,D:1):1);"; append = true)
    summaries = read_trees(path; burnin = 1) do trees
        [evol_distinct(tree).scores for tree in trees]
    end
    sampled = tree_sampler(path; n = 2, replace = false,
        rng = MersenneTwister(42))
    (; summaries, sampled)
end
posterior_example.sampled.sample_indices
posterior_example.sampled.source_indices
```

The returned `TreeStream` is single-pass. Exhaustion and errors close its file;
the do-block also closes it on an early return. Close a manually managed
stream after stopping early. Empty files yield no trees.
Retained trees require nonempty, unique tip labels; parse errors report the
source path, original tree number and line.

`tree_sampler` returns an iterable `TreeSampleResult`. It defaults to `n=100`
with replacement, matching R's sampler; `replace=false` draws distinct records.
Every eligible tree has equal probability, including across files of different
lengths. Pass a seeded `rng` to reproduce Julia draws; R uses a different random
sequence. Empty populations and excessive draws without replacement are errors.

`trees`, `sample_indices`, `source_paths`, `source_indices` and `source_lines`
follow draw order. Sample indices address the filtered population across files;
source indices are original per-file tree numbers. Repeated draws own
independent tree arrays. Sampling scans twice and parses only chosen trees;
unsampled records still undergo framing checks, but not structural parsing.
Keep input files unchanged throughout reading or sampling.

### NEXUS posterior files

A leading `#NEXUS` header selects NEXUS parsing, regardless of the filename
extension. `TREE` and `UTREE` definitions in all `TREES` blocks are read in
file order. `TRANSLATE` applies only to tips in its own block; internal numeric
labels and root labels are preserved. Tree names and bracket annotations,
including `[&R]` and `[&U]`, are discarded; the written topology is not rerooted.

```julia
nexus_example = mktempdir() do dir
    path = joinpath(dir, "posterior.trees")
    write(path, """
    #NEXUS
    BEGIN TREES;
      TRANSLATE 1 A, 2 B, 3 C;
      TREE STATE_0 = [&R] ((1:1,2:2):3,3:4);
      TREE STATE_100 = [&R] ((3:3,1:2):1,2:5);
    END;
    """)
    data = CommunityData(sparse([1 1 0; 1 0 1; 0 1 1]);
        species_labels = ["C", "A", "B"])
    average = read_trees(path) do trees
        mean_dist(data, trees)
    end
    sampled = tree_sampler(path; n = 2, replace = false,
        rng = MersenneTwister(42))
    (; average, sampled, regions = phyloregion(average; k = 2))
end
```

Filtering and `source_indices` count tree definitions across the file's blocks;
metadata commands do not count. `source_lines` points to each `TREE`/`UTREE`
command. Newick and NEXUS paths may be mixed. Tables must precede trees, have
unique keys and taxon labels, and appear at most once per block. Each retained
tip must match a table key; without a table, labels remain literal. Labels are
not renamed, and `TAXA` indices/`LINK` are not resolved (`LINK` is rejected).
Other blocks are skipped, and every block must end with `END` or `ENDBLOCK`.
Command/table validation applies even during burn-in and sampling's counting
scan; tree syntax and tip translation are checked only for parsed trees.

Pass either iterator to the [posterior mean-distance workflow](regionalization.md#mean-distances-across-posterior-trees).

## Prune and Map Trees

```julia
pruned = keep_tip(tree, ["A", "C", "D"])
without_b = drop_tip(tree, "B")

ref = build_phyloref(pruned)
node_map = phyloref_node_table(pruned; map = ref)
edge_map = phyloref_edge_table(pruned; map = ref)

ape_edges = ape_edge_matrix(pruned; map = ref)
signature = phylo_edge_signature(pruned, 1; map = ref)
```

`PhyloRef` provides stable node and edge identities, including conversions
between `CompactTree` IDs and R/ape numbering. Its tables and signatures make
tree results easier to compare across tools and runtimes.

## Taxonomy-aided tree completion

`phylobuilder` explicitly places missing species using original-tree genus
anchors or named internal taxonomic nodes. It returns a `PhyloBuilderResult`
with an independent `tree` and provenance for every requested species.

```julia
using DataFrames

source = mktempdir() do dir
    path = joinpath(dir, "source.nwk")
    write_tree(path, "((G_a:1,G_b:1)G:1,(H_a:1,I_a:1)Family:1)Root;")
    read_tree(path)
end
taxonomy = DataFrame(species = ["G_a", "G_c", "J_a", "Z_a"],
    genus = ["G", "G", "J", "Z"], family = [missing, missing, "Family", missing])
built = phylobuilder(taxonomy, source; extract = false)
built.status             # [:tree, :genus, :family, :missing]
built.added_species      # ["G_c", "J_a"]
built.unmatched_species  # ["Z_a"]
built.tree
```

The species column contains complete exact tip labels, not species epithets.
A label vector is also accepted; genus is inferred only from labels containing
an underscore or space. Optional table columns select genus and higher taxa;
`taxonomy_cols=(:family,:suborder,:order,:class,:phylum)` sets fallback order.
No label spelling is changed automatically.

A missing species attaches to the original congeneric MRCA, or its parent
when only one congener exists. If no congeneric tip exists, try a named
internal genus node, then higher taxonomic nodes. All anchors come from the
original tree; additions cannot seed later placements. Nonmonophyletic genera
fail unless `nonmonophyly=:mrca` is explicit. Repeated matching internal names
are ambiguous and fail. Unmatched species remain reported with `:missing`.

By default new tip lengths use the anchor's longest descendant path and the
source must be ultrametric (`atol=0`, `rtol=1e-8`). For other trees explicitly
set `branch_length=:zero` or a finite nonnegative terminal length. Existing
edges retain their lengths. Default `extract=true` then prunes to the matched
requested species, which can contract unary paths and remove node labels;
`false` retains every original tip. `anchor_nodes` records original compact
node IDs even though output IDs/caches are rebuilt. `align_community` still
rejects unknown species; apply this completion step explicitly when intended.

## Tree branch time slices

`timeslice(tree, depth; collapse=false)` removes the most recent `depth`
units of branch history from a dated ultrametric tree. It returns a new
`CompactTree`, retaining all species in their original column order.
Older branches retain their lengths, crossing branches are shortened, and
branches starting at or younger than the cut become zero, including at exact
node boundaries. Depth must be between zero and the tree height, inclusive;
zero preserves lengths and the full height makes every edge zero.

```julia
mktempdir() do dir
    path = joinpath(dir, "tree.nwk")
    write_tree(path, "(((A:1,B:1):1,C:2):2,(D:2,E:2):2);")
    tree = read_tree(path)
    data = CommunityData(sparse([1 1 0 0 0; 0 0 1 0 0; 0 0 0 1 1]);
        species_labels = tree.tip_labels)
    profiles = map((0.0, 1.0, 2.0)) do depth
        sliced = timeslice(tree, depth; collapse = true)
        ctx = prepare_community(data, sliced)
        (; depth, pe = phylo_endemism(ctx).pe,
            regions = phyloregion(phylobeta(ctx); k = 2))
    end
end
```

By default topology, node IDs, edge order and labels are retained, with root
distances rebuilt. `collapse=true` contracts exactly zero internal edges,
including pre-existing ones, while retaining the root and every tip. IDs are
then rebuilt and labels on removed nodes disappear. Positive edges are kept
regardless of size. Rooted PD/PE and length-based beta agree before and after
collapse; edge-indexed diagnostics and topology-based statistics may change.

Validation uses the same `atol=0.0`, `rtol=1e-8` and longest descendant ages as
`get_clades`. For accepted near-ultrametric trees, cuts use these local parent
ages; tolerances do not move the cut or collapse small positive edges. The
operation takes O(nodes+tips) time/storage and never changes its input.
Use the original tree for each absolute depth and prepare a new community
context per sliced tree. At full height all PD/PE values are zero and beta is
undefined (`NaN`), so that endpoint cannot be passed directly to regionalization.
