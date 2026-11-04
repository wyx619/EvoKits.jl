# Per-file NEXUS state. Translation tables belong to one TREES block and are
# never retained per posterior draw.
mutable struct _NexusTreeState
    block::String
    block_line::Int
    translation::Dict{String, String}
    translated::Bool
    tree_seen::Bool
end

_NexusTreeState() = _NexusTreeState("", 0, Dict{String, String}(), false, false)

function _nexus_reset!(state)
    state.block, state.block_line = "", 0
    state.translation = Dict{String, String}()
    state.translated, state.tree_seen = false, false
    return nothing
end

const _NEXUS_DELIMITERS = ('(', ')', ',', ':', ';', '=', '*')

# Only command headers and TRANSLATE tables are tokenized. The Newick body is
# left intact, so skipped posterior trees do not undergo structural parsing.
function _nexus_token(text, index)
    while index <= lastindex(text) && isspace(text[index])
        index = nextind(text, index)
    end
    index > lastindex(text) && return nothing
    start, c = index, text[index]
    if c == '\''
        buffer = IOBuffer()
        index = nextind(text, index)
        while index <= lastindex(text)
            c = text[index]
            index = nextind(text, index)
            if c == '\''
                if index <= lastindex(text) && text[index] == '\''
                    write(buffer, '\'')
                    index = nextind(text, index)
                else
                    index <= lastindex(text) && !isspace(text[index]) &&
                        !(text[index] in _NEXUS_DELIMITERS) &&
                        throw(ArgumentError("unexpected text after a quoted NEXUS token"))
                    return (value=String(take!(buffer)), next=index, quoted=true)
                end
            else
                write(buffer, c)
            end
        end
        throw(ArgumentError("unterminated quoted NEXUS token"))
    end
    if c in _NEXUS_DELIMITERS
        return (value=string(c), next=nextind(text, index), quoted=false)
    end
    while index <= lastindex(text) && !isspace(text[index]) && !(text[index] in _NEXUS_DELIMITERS)
        text[index] in ('\'', '"') && throw(ArgumentError("use single quotes around complete NEXUS labels"))
        index = nextind(text, index)
    end
    return (value=String(SubString(text, start, prevind(text, index))), next=index, quoted=false)
end

_nexus_is(token, value) = token !== nothing && !token.quoted && token.value == value

function _nexus_label(token)
    token !== nothing && !isempty(strip(token.value)) &&
        (token.quoted || !(first(token.value) in _NEXUS_DELIMITERS)) ||
        throw(ArgumentError("expected a nonempty NEXUS name"))
    return token.value
end

function _nexus_end(text, index)
    token = _nexus_token(text, index)
    _nexus_is(token, ";") && _nexus_token(text, token.next) === nothing ||
        throw(ArgumentError("unexpected text at the end of a NEXUS command"))
    return nothing
end

function _nexus_translate!(state, text, index)
    state.translated && throw(ArgumentError("only one TRANSLATE table is allowed per TREES block"))
    state.tree_seen && throw(ArgumentError("TRANSLATE must precede every tree in its TREES block"))
    labels = Set{String}()
    while true
        key = _nexus_token(text, index)
        name = _nexus_label(key)
        value = _nexus_token(text, key.next)
        label = _nexus_label(value)
        haskey(state.translation, name) && throw(ArgumentError("duplicate TRANSLATE key: $name"))
        label in labels && throw(ArgumentError("duplicate TRANSLATE taxon label: $label"))
        state.translation[name] = label
        push!(labels, label)
        separator = _nexus_token(text, value.next)
        if _nexus_is(separator, ";")
            _nexus_end(text, value.next)
            break
        end
        _nexus_is(separator, ",") || throw(ArgumentError("expected a comma between TRANSLATE entries"))
        index = separator.next
    end
    state.translated = true
    return nothing
end

function _nexus_tree_record!(stream, record)
    state, text = stream.nexus, record.text
    token = _nexus_token(text, firstindex(text))
    token === nothing && return nothing
    token.quoted && throw(ArgumentError("NEXUS command names must not be quoted"))
    command = uppercase(token.value)
    if command == "BEGIN"
        isempty(state.block) || throw(ArgumentError("nested NEXUS blocks are unsupported"))
        name = _nexus_token(text, token.next)
        block = lowercase(_nexus_label(name))
        _nexus_end(text, name.next)
        state.block, state.block_line = block, record.line
        return nothing
    elseif command in ("END", "ENDBLOCK")
        isempty(state.block) && throw(ArgumentError("END outside a NEXUS block"))
        _nexus_end(text, token.next)
        _nexus_reset!(state)
        return nothing
    end
    isempty(state.block) && throw(ArgumentError("expected BEGIN before a NEXUS command"))
    state.block == "trees" || return nothing
    if command == "TRANSLATE"
        _nexus_translate!(state, text, token.next)
        return nothing
    elseif command == "TITLE"
        title = _nexus_token(text, token.next)
        _nexus_label(title)
        _nexus_end(text, title.next)
        return nothing
    elseif !(command in ("TREE", "UTREE"))
        throw(ArgumentError("unsupported TREES command: $command (supported: TREE, UTREE, TRANSLATE, TITLE)"))
    end
    name = _nexus_token(text, token.next)
    _nexus_is(name, "*") && (name = _nexus_token(text, name.next))
    _nexus_label(name)
    separator = _nexus_token(text, name.next)
    _nexus_is(separator, "=") || throw(ArgumentError("expected '=' after the tree name"))
    index = separator.next
    while index <= lastindex(text) && isspace(text[index])
        index = nextind(text, index)
    end
    index <= lastindex(text) && text[index] != ';' || throw(ArgumentError("empty TREE definition"))
    state.tree_seen = true
    return (text=String(SubString(text, index)), line=record.line,
        translation=state.translated ? state.translation : nothing)
end
