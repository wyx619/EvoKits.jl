"""
    TreeStream

Single-pass, closeable iterator of `CompactTree`s created by [`read_trees`](@ref).
It owns at most one open file and one Newick record at a time. Use the do-block
form of `read_trees` when iteration may stop early.
"""
mutable struct TreeStream
    paths::Vector{String}
    burnin::Int
    thin::Int
    file_index::Int
    io::Union{Nothing, IOStream}
    record_index::Int
    line::Int
    closed::Bool
end

Base.IteratorSize(::Type{TreeStream}) = Base.SizeUnknown()
Base.eltype(::Type{TreeStream}) = CompactTree
Base.isopen(stream::TreeStream) = !stream.closed

function Base.close(stream::TreeStream)
    stream.io === nothing || close(stream.io)
    stream.io = nothing
    stream.closed = true
    return nothing
end

function _tree_stream_integer(value, name, minimum)
    value isa Integer && !(value isa Bool) && minimum <= value <= typemax(Int) ||
        throw(ArgumentError("$name must be an integer >= $minimum representable as Int"))
    return Int(value)
end

function _tree_stream_paths(input)
    paths = input isa AbstractString ? [input] : input
    paths isa Union{AbstractVector, Tuple} && !isempty(paths) &&
        all(path -> path isa AbstractString, paths) ||
        throw(ArgumentError("supply a Newick file path or a nonempty vector/tuple of file paths"))
    result = String[abspath(path) for path in paths]
    all(isfile, result) || throw(ArgumentError("every tree path must name an existing file"))
    identities = Sys.iswindows() ? lowercase.(realpath.(result)) : realpath.(result)
    allunique(identities) || throw(ArgumentError("tree file paths must not repeat"))
    return result
end

"""
    read_trees(paths; burnin=0, thin=1)
    read_trees(f, paths; burnin=0, thin=1)

Read semicolon-terminated Newick records one at a time, returning a `TreeStream`
of `CompactTree`s. Supply a file path or a nonempty vector/tuple of paths;
files are concatenated in the supplied order. Repeated file paths are errors.
Records may share a line or span lines. Whitespace, an initial UTF-8 BOM,
square-bracket comments (including rooting annotations), and single-quoted
labels with doubled quote escapes are supported. NEXUS containers are rejected.

In each file, discard the first `burnin` records, then retain records
`burnin+1`, `burnin+1+thin`, ... . Empty files yield no trees. Framing errors
such as unterminated records, quotes or comments fail, also in skipped records;
tree syntax is parsed only for retained records. Errors include file, original
one-based tree number and line. Missing branch lengths remain undefined as in
`read_tree`; downstream metrics apply their own branch-length requirements.

The stream is single-pass and opens files lazily. Exhaustion and read/parse
errors close the current file. Call `close(stream)` after an early stop, or
use `read_trees(paths) do trees ... end` for deterministic cleanup. File
contents must remain unchanged during iteration. No collection of all tree
strings or parsed trees is retained.
"""
function read_trees(paths; burnin=0, thin=1)
    files = _tree_stream_paths(paths)
    discarded = _tree_stream_integer(burnin, "burnin", 0)
    stride = _tree_stream_integer(thin, "thin", 1)
    stream = TreeStream(files, discarded, stride, 0, nothing, 0, 1, false)
    finalizer(close, stream)
    return stream
end

function read_trees(f::Function, paths; kwargs...)
    stream = read_trees(paths; kwargs...)
    try
        return f(stream)
    finally
        close(stream)
    end
end

function _tree_stream_error(stream, line, message)
    path = stream.paths[stream.file_index]
    return ArgumentError("Newick file $path, tree $(stream.record_index+1), line $line: $message")
end

function _tree_stream_open_next!(stream)
    stream.io === nothing || close(stream.io)
    stream.io = nothing
    stream.file_index += 1
    if stream.file_index > length(stream.paths)
        close(stream)
        return false
    end
    stream.io = open(stream.paths[stream.file_index], "r")
    stream.record_index, stream.line = 0, 1
    return true
end

# Record framing respects quotes and nested comments, so metadata or labels
# containing semicolons cannot silently split one tree into several records.
function _tree_stream_record!(stream)
    buffer = IOBuffer()
    quoted, comments, started = false, 0, false
    start_line = stream.line
    while !eof(stream.io)
        at_start = position(stream.io) == 0
        c = read(stream.io, Char)
        line = stream.line
        c == '\n' && (stream.line += 1)
        at_start && c == '\ufeff' && continue
        if comments > 0
            c == '[' && (comments += 1)
            c == ']' && (comments -= 1)
            continue
        end
        if !quoted && c == '['
            comments = 1
            started || (start_line = line)
            continue
        end
        !quoted && c == ']' && throw(_tree_stream_error(stream, line, "unmatched closing comment bracket"))
        if !started
            isspace(c) && continue
            started, start_line = true, line
        end
        write(buffer, c)
        c == '\'' && (quoted = !quoted)
        !quoted && c == ';' && return (text=String(take!(buffer)), line=start_line)
    end
    comments == 0 || throw(_tree_stream_error(stream, start_line, "unterminated square-bracket comment"))
    !quoted || throw(_tree_stream_error(stream, start_line, "unterminated quoted label"))
    started && throw(_tree_stream_error(stream, start_line, "record is missing its trailing semicolon"))
    return nothing
end

function _next_tree_record!(stream)
    stream.closed && return nothing
    while true
        stream.io === nothing && !_tree_stream_open_next!(stream) && return nothing
        record = _tree_stream_record!(stream)
        if record === nothing
            _tree_stream_open_next!(stream) || return nothing
            continue
        end
        stream.record_index == 0 && occursin(r"^#NEXUS"i, record.text) && throw(_tree_stream_error(
            stream, record.line, "NEXUS input is unsupported; export trees as Newick"))
        stream.record_index += 1
        index = stream.record_index
        if index > stream.burnin && (index - stream.burnin - 1) % stream.thin == 0
            return (; record.text, path=stream.paths[stream.file_index], index, record.line)
        end
    end
end

# Validate Newick structure before invoking the permissive dependency parser.
# Replace quoted labels by collision-free plain tokens and restore them after
# parsing; the dependency itself does not decode quoted Newick labels.
function _tree_stream_normalize(text)
    prefix = "__EVOKITS_QUOTED_"
    while occursin(prefix, text)
        prefix *= "_"
    end
    labels = Dict{String, String}()
    output, token = IOBuffer(), IOBuffer()
    index, depth, state = firstindex(text), 0, :subtree
    delimiters = ('(', ')', ',', ':', ';')
    while index <= lastindex(text)
        c = text[index]
        if isspace(c)
            index = nextind(text, index)
            continue
        end
        if c in delimiters
            if c == '('
                state === :subtree || throw(ArgumentError("unexpected opening parenthesis"))
                depth += 1
            elseif c == ')'
                state in (:suffix, :internal, :after_length) && depth > 0 ||
                    throw(ArgumentError("unexpected closing parenthesis or missing child"))
                depth -= 1
                state = :internal
            elseif c == ','
                state in (:suffix, :internal, :after_length) && depth > 0 ||
                    throw(ArgumentError("unexpected comma or missing child"))
                state = :subtree
            elseif c == ':'
                state in (:suffix, :internal) || throw(ArgumentError("unexpected branch-length separator"))
                state = :length
            else
                state in (:suffix, :internal, :after_length) && depth == 0 ||
                    throw(ArgumentError("incomplete tree at semicolon"))
                state = :done
            end
            write(output, c)
            index = nextind(text, index)
            continue
        end
        quoted = c == '\''
        if quoted
            index = nextind(text, index)
            closed = false
            while index <= lastindex(text)
                c = text[index]
                index = nextind(text, index)
                if c == '\''
                    if index <= lastindex(text) && text[index] == '\''
                        write(token, '\'')
                        index = nextind(text, index)
                    else
                        closed = true
                        break
                    end
                else
                    write(token, c)
                end
            end
            closed || throw(ArgumentError("unterminated quoted label"))
        else
            while index <= lastindex(text) && !(text[index] in delimiters)
                text[index] == '\'' && throw(ArgumentError("quotes must enclose the complete label"))
                write(token, text[index])
                index = nextind(text, index)
            end
        end
        value = String(take!(token))
        quoted || (value = strip(value))
        if state === :length
            !quoted && tryparse(Float64, value) !== nothing || throw(ArgumentError("invalid branch length"))
            state = :after_length
        elseif state in (:subtree, :internal)
            isempty(value) && throw(ArgumentError("empty node label"))
            state = :suffix
        else
            throw(ArgumentError("unexpected text after node label or branch length"))
        end
        if quoted
            placeholder = prefix * string(length(labels) + 1)
            labels[placeholder] = value
            value = placeholder
        end
        write(output, value)
    end
    state === :done || throw(ArgumentError("incomplete Newick record"))
    return String(take!(output)), labels
end

function _parse_tree_record(record)
    try
        newick, labels = _tree_stream_normalize(record.text)
        node = if startswith(newick, '(')
            NewickTree.readnw(newick, UInt32)
        else
            # The dependency parser leaves a standalone leaf's data unset.
            label = first(split(newick, [':',';']; limit=2))
            NewickTree.Node(UInt32(1), NewickData(d=0.0, n=String(label)))
        end
        tree = compact_tree(node)
        if !isempty(labels)
            for names in (tree.node_labels, tree.tip_labels), i in eachindex(names)
                names[i] = get(labels, names[i], names[i])
            end
            empty!(tree.tipname_to_id)
            for i in eachindex(tree.tip_labels)
                tree.tipname_to_id[tree.tip_labels[i]] = tree.tip_ids[i]
            end
        end
        _validated_labels(tree.tip_labels, tree.ntips, "tree tip")
        return tree
    catch error
        error isa InterruptException && rethrow()
        message = first(sprint(showerror, error), 240)
        throw(ArgumentError("Newick file $(record.path), tree $(record.index), line $(record.line): $message"))
    end
end

function Base.iterate(stream::TreeStream, state=nothing)
    try
        record = _next_tree_record!(stream)
        record === nothing && return nothing
        return _parse_tree_record(record), nothing
    catch
        close(stream)
        rethrow()
    end
end
