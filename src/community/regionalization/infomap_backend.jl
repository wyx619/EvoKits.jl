# IGraphs 1.x wraps the igraph 1.x ABI. Keep native objects and RNG switching
# inside this adapter; no Graphs graph or dense adjacency matrix is needed.
const _INFOMAP_BACKEND_LOCK = ReentrantLock()

function _infomap_native_check(code, operation)
    code == IGraphs.LibIGraph.IGRAPH_SUCCESS ||
        error("igraph $operation failed with error code $(Int(code))")
    return nothing
end

function _infomap_engine_version()
    text = Ref{Ptr{Cchar}}()
    major, minor, patch = Ref{Cint}(), Ref{Cint}(), Ref{Cint}()
    IGraphs.LibIGraph.igraph_version(text, major, minor, patch)
    return unsafe_string(text[])
end

function _infomap_backend(edges, edge_weights, nvertices, trials, seed,
    regularized, regularization_strength)
    # IGraphs' error handler has shared state. Serialize our native calls and
    # their finalizers, including graph construction and destruction.
    return lock(_INFOMAP_BACKEND_LOCK) do
        lib = IGraphs.LibIGraph
        graph = IGraphs.IGraph(nvertices)
        native_edges = IGraphs.IGVectorInt()
        native_weights = IGraphs.IGVectorFloat()
        native_membership = IGraphs.IGVectorInt()
        try
            lib.vector_int_resize(native_edges, length(edges))
            copyto!(native_edges, edges)
            _infomap_native_check(lib.igraph_add_edges(graph.objref,
                native_edges.objref, C_NULL), "edge construction")
            if edge_weights !== nothing
                lib.vector_resize(native_weights, length(edge_weights))
                copyto!(native_weights, edge_weights)
            end

            # A private PCG32 state makes the seed independent of both Julia's
            # RNG and the RNG selected by other users of igraph. igraph 1.x
            # set_default returns the previous pointer, which must be restored
            # before the private state is destroyed (including on exceptions).
            native_rng = Ref{IGraphs.LibIGraph.igraph_rng_t}()
            rng_type = cglobal((:igraph_rngtype_pcg32,
                IGraphs.LibIGraph.libigraph), IGraphs.LibIGraph.igraph_rng_type_t)
            _infomap_native_check(lib.igraph_rng_init(native_rng, rng_type), "RNG initialization")
            try
                _infomap_native_check(lib.igraph_rng_seed(native_rng, seed), "RNG seeding")
                GC.@preserve native_rng graph native_edges native_weights native_membership begin
                    previous_rng = lib.igraph_rng_set_default(native_rng)
                    codelength = try
                        only(lib.community_infomap(graph,
                            edge_weights === nothing ? IGraphs.IGNull() : native_weights,
                            IGraphs.IGNull(), trials, regularized,
                            regularization_strength, native_membership))
                    finally
                        lib.igraph_rng_set_default(previous_rng)
                    end
                    membership = Vector(native_membership)
                    length(membership) == nvertices && all(>=(0), membership) ||
                        error("igraph returned an invalid Infomap membership")
                    isfinite(codelength) && codelength >= 0 ||
                        error("igraph returned a nonfinite or negative Infomap codelength; rescale weights or regularization strength")
                    return membership, codelength, _infomap_engine_version()
                end
            finally
                lib.igraph_rng_destroy(native_rng)
            end
        finally
            # Release native allocations promptly instead of waiting for GC.
            # These are the Ref objects on which IGraphs installs finalizers.
            finalize(native_membership.objref)
            finalize(native_weights.objref)
            finalize(native_edges.objref)
            finalize(graph.objref)
        end
    end
end
