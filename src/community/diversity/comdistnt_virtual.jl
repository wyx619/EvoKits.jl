# Each unique unordered pair needs both directional means (one for a self
# pair). Group those queries by target to reuse each nearest-target pass.
function _comdistnt_pair_plan(rows, selected)
    mapping = zeros(Int, length(selected))
    jobs = NTuple{3, Int}[]
    used = falses(length(rows.nodes))
    previous, npairs = (0, 0), 0
    for output in sortperm(selected; by = pair -> minmax(pair...))
        i, j = minmax(selected[output]...)
        (isempty(rows.nodes[i]) || isempty(rows.nodes[j])) && continue
        if (i, j) != previous
            npairs += 1
            push!(jobs, (i, j, npairs))
            i != j && push!(jobs, (j, i, npairs))
            used[i] = used[j] = true
            previous = (i, j)
        end
        mapping[output] = npairs
    end
    sort!(jobs)
    return (; mapping, jobs, used, npairs)
end

function _comdistnt_selected(rows, tree, plan, weighted, exclude)
    values = fill(NaN, plan.npairs)
    workspace = _ComdistntWorkspace(tree.nnodes)
    previous = 0
    for (target, source, pair) in plan.jobs
        if target != previous
            _comdistnt_target!(workspace, tree, rows.nodes[target], rows.lengths)
            previous = target
        end
        value = _comdistnt_directional(workspace, rows.nodes[source], rows.raw[source], exclude)
        # The smaller target is visited first. Store source j -> target i,
        # then combine i -> j using the original source-total weighting.
        values[pair] = target <= source ? value :
            _comdistnt_symmetric(value, values[pair], rows, source, target, weighted)
    end
    return [pair == 0 ? NaN : values[pair] for pair in plan.mapping]
end

function _comdistnt_adaptive(rows, tree, selected, weighted, exclude; strategy = :auto)
    strategy in (:auto, :full, :virtual) || throw(ArgumentError("unknown private comdistnt strategy"))
    if selected === nothing
        targets = count(!isempty, rows.nodes)
        targets == 0 && return fill(NaN, length(rows.nodes), length(rows.nodes))
        plan = nothing
    else
        isempty(selected) && return Float64[]
        plan = _comdistnt_pair_plan(rows, selected)
        targets = count(plan.used)
        targets == 0 && return fill(NaN, length(selected))
        # Only endpoints belonging to a nonempty pair affect any result.
        if any(i -> !plan.used[i] && !isempty(rows.nodes[i]), eachindex(plan.used))
            rows = merge(rows, (; nodes = [plan.used[i] ? row : Int32[] for (i, row) in enumerate(rows.nodes)]))
        end
    end
    virtual = strategy === :virtual || (strategy === :auto && _comdist_use_virtual(rows, tree, targets))
    if virtual
        rows, tree = _comdist_virtual(rows, tree)
    elseif isempty(rows.lengths)
        rows = merge(rows, (; lengths = _comdist_lengths(tree)))
    end
    return selected === nothing ? _comdistnt(rows, tree, weighted, exclude) :
        _comdistnt_selected(rows, tree, plan, weighted, exclude)
end
